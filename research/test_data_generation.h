#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>

using std::vector;
using std::pair;
using std::string;
using std::time_t;
using std::transform;
using std::uniform_real_distribution;
using std::mt19937;
using std::random_device;
using std::back_inserter;
using std::lower_bound;
using std::distance;

// Типы синонимов для улучшения читаемости кода
using TimeVector = vector<time_t>;
using ParamVector = vector<double>;
using ParamPair = pair<TimeVector, ParamVector>;
//SynteticTimeSeriesGenerator
class SynteticTimeSeriesGenerator {
public:
    SynteticTimeSeriesGenerator(const vector<pair<string, double>>& params, double default_duration = 350000,
        pair<double, double> timeDist = { 200, 400 },
        pair<double, double> valueDist = { 0.9998, 1.0002 })
        : parameters(params), Duration(default_duration), timeDistribution(timeDist), valueDistribution(valueDist), gen(rd()) {
        uniform_real_distribution<double> timeDis(timeDistribution.first, timeDistribution.second);

        for (const auto& param : parameters) {
            TimeVector timeValues;
            ParamVector paramValues;

            uniform_real_distribution<double> normalDis(param.second * valueDistribution.first, param.second * valueDistribution.second);

            double timeStep = timeDis(gen);
            for (double time = std::time(nullptr); time <= std::time(nullptr) + Duration; time += timeStep) {
                timeValues.push_back(static_cast<time_t>(time));
                timeStep = timeDis(gen);
            }

            transform(timeValues.begin(), timeValues.end(), back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }

    void applyJump(double jumpTime, double jumpValue, const string& paramName) {
        for (size_t i = 0; i < parameters.size(); ++i) {
            if (parameters[i].first == paramName) {
                auto it = lower_bound(data[i].first.begin(), data[i].first.end(), time(nullptr) + jumpTime);
                size_t position = distance(data[i].first.begin(), it);
                uniform_real_distribution<double> normalDis(jumpValue * valueDistribution.first, jumpValue * valueDistribution.second);
                for (size_t j = position; j < data[i].first.size(); ++j) {
                    double value = normalDis(gen);
                    data[i].second[j] = value;
                }
            }
        }
    }

    vector<ParamPair> getData() const {
        return data;
    }

private:
    vector<pair<string, double>> parameters;
    vector<ParamPair> data;
    random_device rd;
    mt19937 gen;
    pair<double, double> timeDistribution;
    pair<double, double> valueDistribution;
    double Duration;
};

TEST(Random, PrepareTimeSeries)
{
    const vector<pair<string, double>> parameters = {
        { "rho_in", 860 },
        { "visc_in", 15e-6},
        { "p_in", 6e6},
        { "Q", 0.2 },
    };

    SynteticTimeSeriesGenerator dataGenerator(parameters);
    

    const double jumpTime = 100000;
    const double jumpValue = 870;
    dataGenerator.applyJump(jumpTime, jumpValue, "rho_in");

    const auto data = dataGenerator.getData();

    vector_timeseries_t params(data);

    // Задаём интересующий нас момент времени
    time_t test_time = StringToUnix("24.03.2024 08:53:50");

    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
}

#pragma once

/// @brief Тесты для солвера quickest_ultimate_fv_solver
class QuickWithQuasiStationaryModel : public ::testing::Test {
protected:

protected:
    /// @brief Параметры трубы
    pipe_properties_t pipe;
    /// @brief Профиль расхода
    vector<double> Q_profile;
    /// @brief Модель адвекции
    std::unique_ptr<PipeQAdvection> advection_model;
protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 200км, с шагом разбиения для расчтной сетки 100 м, диаметром 514мм
        simple_pipe_properties simple_pipe;
        simple_pipe.length = 200e3;
        simple_pipe.diameter = 0.514;
        simple_pipe.dx = 100;

        pipe = pipe_properties_t::build_simple_pipe(simple_pipe);
        Q_profile = vector<double>(pipe.profile.getPointCount(), 0.2); // задаем по трубе расход 0.2 м3/с
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
    }
};

/// @brief Слой для расчета плотности, вязкости методом конечных объемов 
struct density_viscosity_cell_layer {
    /// @brief Профиль плотности
    std::vector<double> density;
    /// @brief Профиль вязкости
    std::vector<double> viscosity;
    /// @brief Профиль вспомогательных расчетов для метода конечных объемов (и для вязкости, и для плотности)
    quickest_ultimate_fv_solver_traits<1>::specific_layer specific;
    /// @brief Инициализация профилей
    /// @param point_count Количество точек
    density_viscosity_cell_layer(size_t point_count)
        : density(point_count - 1)
        , viscosity(point_count - 1)
        , specific(point_count)
    {
    }

    // @brief Подготовка плотности для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_density_quick_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.density, layer.specific);
    }
    /// @brief Подготовка вязкости для расчета методом конечных объемов 
    static quickest_ultimate_fv_wrapper<1> get_viscosity_quick_wrapper(density_viscosity_cell_layer& layer)
    {
        return quickest_ultimate_fv_wrapper<1>(layer.viscosity, layer.specific);
    }
};

/// @brief Уравнение трубы для задачи PQ
class pipe_model_PQ_cell_parties_t : public ode_t<1>
{
public:
    using ode_t<1>::equation_coeffs_type;
    using ode_t<1>::right_party_type;
    using ode_t<1>::var_type;
protected:
    const vector<double>& rho_profile;
    const vector<double>& nu_profile;
    const pipe_properties_t& pipe;
    const double flow;
    const int solver_direction;
public:
    /// @brief Констуктор уравнения трубы
    /// @param pipe Ссылка на сущность трубы
    /// @param oil Ссылка на сущность нефти
    /// @param flow Объемный расход
    /// @param solver_direction Направление расчета по Эйлеру, должно обязательно совпадать с параметром солвера Эйлера
    pipe_model_PQ_cell_parties_t(pipe_properties_t& pipe, vector<double>& rho_profile, vector<double>& nu_profile, double flow,
        int solver_direction)
        : pipe(pipe)
        , rho_profile(rho_profile)
        , nu_profile(nu_profile)
        , flow(flow)
        , solver_direction(solver_direction)
    {

    }

    /// @brief Возвращает известную уравнению сетку
    virtual const vector<double>& get_grid() const override {
        return pipe.profile.coordinates;
    }

    /// @brief Возвращает значение правой части ДУ
    /// @param grid_index Обсчитываемый индекс расчетной сетки
    /// @param point_vector Начальные условия
    /// @return Значение правой части ДУ в точке point_vector
    virtual right_party_type ode_right_party(
        size_t grid_index, const var_type& point_vector) const override
    {

        /// Обработка индекса в случае расчетов на границах трубы
        /// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
        //grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
        //grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

        size_t reo_index = solver_direction == +1
            ? grid_index
            : grid_index - 1;

        double rho = rho_profile[reo_index];
        double S_0 = pipe.wall.getArea();
        double v = flow / (S_0);
        double Re = v * pipe.wall.diameter / nu_profile[reo_index];
        double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
        double tau_w = lambda / 8 * rho * v * abs(v);

        double height_derivative = pipe.profile.get_height_derivative(grid_index, solver_direction);
        double result = -4 * tau_w / pipe.wall.diameter - rho * M_G * height_derivative;
        return result;
    }
};

TEST_F(QuickWithQuasiStationaryModel, WorkingWithTimeSeries)
{
    const vector<pair<string, double>> parameters = {
        { "rho_in", 860 },
        { "visc_in", 15e-6},
        { "p_in", 6e6},
        { "Q", 0.2 },
    };

    const double duration = 350000;

    SynteticTimeSeriesGenerator dataGenerator(parameters);

    const auto data = dataGenerator.getData();

    //vector_timeseries_t params(data);

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0];
    double v = advection_model->getEquationsCoeffs(0, 0);

    ring_buffer_t<density_viscosity_cell_layer> buffer(2, pipe.profile.getPointCount());


    //double v_max = пересчитываю максимальное значение из временного ряда расхода на скорость
    double dt = abs(dx);// double dt = abs(dx/v_max);

    auto& rho_initial = buffer[0].density;
    auto& viscosity_initial = buffer[0].viscosity;
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // инициализация начальной плотности

    buffer.advance(+1);

    double t = std::time(nullptr);

    do
    {

    } while (t < t + duration - dt);
}