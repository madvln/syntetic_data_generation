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

/// @brief Функция расчета скорости из расхода
/// @param Q Расход, (м^3/с)
/// @param internal_diameter Внутренний диаметр, (м)
/// @return Скорость, (м/с)
double calc_speed(double Q, double internal_diameter)
{
    double speed = (4 * Q) / (pow(internal_diameter, 2) * M_PI);
    return speed;
}

// Типы синонимов для улучшения читаемости кода
using TimeVector = vector<time_t>;
using ParamVector = vector<double>;
using ParamPair = pair<TimeVector, ParamVector>;


struct timeseries_generator_settings {
    std::time_t duration{ 350000 };
    std::time_t sample_time_min{200};
    std::time_t sample_time_max{400};
    double value_relative_increment{ 0.0002 };
    double value_relative_decrement{ 0.0002 };
    vector<pair<string, double>> timeseries_initial_values;
};

/// @brief сделать snake case
class SynteticTimeSeriesGenerator {
public:
    SynteticTimeSeriesGenerator(const vector<pair<string, double>>& params, double default_duration = 350000,
        pair<double, double> timeDist = { 200, 400 },
        pair<double, double> valueDist = { 0.9998, 1.0002 })
        : parameters(params)
        , Duration(default_duration)
        , timeDistribution(timeDist)
        , valueDistribution(valueDist)
        , gen(rd()) 
    {
        std::uniform_real_distribution<double> timeDis(timeDistribution.first, timeDistribution.second);

        for (const auto& param : parameters) {
            TimeVector timeValues;
            ParamVector paramValues;

            std::uniform_real_distribution<double> normalDis(param.second * valueDistribution.first, param.second * valueDistribution.second);

            double timeStep = timeDis(gen);
            for (double time = std::time(nullptr); time <= std::time(nullptr) + Duration; time += timeStep) {
                timeValues.push_back(static_cast<time_t>(time));
                timeStep = timeDis(gen);
            }

            std::transform(timeValues.begin(), timeValues.end(), std::back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }

    void applyJump(double jumpTime, double jumpValue, const string& paramName) {
        for (size_t i = 0; i < parameters.size(); ++i) {
            if (parameters[i].first == paramName) {
                auto it = std::lower_bound(data[i].first.begin(), data[i].first.end(), time(nullptr) + jumpTime);
                size_t position = std::distance(data[i].first.begin(), it);
                std::uniform_real_distribution<double> normalDis(jumpValue * valueDistribution.first, jumpValue * valueDistribution.second);
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
    std::random_device rd;
    std::mt19937 gen;
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

    const double duration = 350000;

    // Задаём настроечные параметры
    pair<double, double> timeDistribution = { 200, 400 };//Размах во времени, от 200 с до 400 с
    pair<double, double> valueDistribution = { 0.9998, 1.0002 };//Процентное отклонение, имитация шума у датчиков
    
    SynteticTimeSeriesGenerator dataGenerator(parameters, duration, timeDistribution, valueDistribution);

    const double jumpTime_rho = 100000;
    const double jumpValue_rho = 870;
    dataGenerator.applyJump(jumpTime_rho, jumpValue_rho, "rho_in");

    const double jumpTime_Q = 150000;
    const double jumpValue_Q = 0.1;
    dataGenerator.applyJump(jumpTime_Q, jumpValue_Q, "Q");

    const auto data = dataGenerator.getData();

    vector_timeseries_t params(data);

    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);

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
    {}

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
    // Для генерации временных рядов
    const vector<pair<string, double>> parameters = {
        { "rho_in", 860 },
        { "visc_in", 15e-6},
        { "p_in", 6e6},
        { "Q", 0.2 },
    };

    const double duration = 300000;// Задаю время моделирования (опционально, по умолчанию 350000)
    
    SynteticTimeSeriesGenerator dataTimeSeries(parameters, duration); // Генерируем данные

    const double jumpTime_Q = 50000; // Момент скачка по расходу
    const double jumpValue_Q = 0.1; // Значение расхода в момент скачка
    dataTimeSeries.applyJump(jumpTime_Q, jumpValue_Q, "Q");
    
    // Получаем данные
    const auto data = dataTimeSeries.getData(); 
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(data);

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0]; // Шаг сетки

    ring_buffer_t<density_viscosity_cell_layer> buffer(2, pipe.profile.getPointCount());

    auto& rho_initial = buffer[0].density;
    auto& viscosity_initial = buffer[0].viscosity;
    rho_initial = vector<double>(rho_initial.size(), 850); // Инициализация начальной плотности
    viscosity_initial = vector<double>(viscosity_initial.size(), 1e-5); // Инициализация начальной вязкости
    double p_initial = 6e6; // Давление на входе изначальное
    buffer.advance(+1);

    double v_max = 1; // Предполагаем скорость для Куранта = 1, скорость, больше чем во временных рядах и в профиле
    double dt = abs(dx/v_max); // Постоянный шаг по времени для Куранта = 1

    double t = std::time(nullptr); // Время на данный момент

    vector<double> initial_p_profile;
    vector<double> diff_p_profile = vector<double>(pipe.profile.getPointCount(), 0);

    do
    {
        vector<double> p_profile(pipe.profile.getPointCount());

        auto density_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_density_quick_wrapper);

        auto viscosity_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_viscosity_quick_wrapper);

        int euler_direction = +1;
        
        if (t == std::time(nullptr))
        {
            pipe_model_PQ_cell_parties_t pipeModel(pipe, density_wrapper.previous().vars, viscosity_wrapper.previous().vars, Q_profile[0], euler_direction);
            solve_euler<1>(pipeModel, euler_direction, p_initial, &p_profile);
            initial_p_profile = p_profile;// получаем изначальный профиль
        }
        t += dt;
        
        // Задаём интересующий нас момент времени
        time_t time_model = static_cast<time_t>(t);
        // Интерополируем значения параметров в заданный момент времени
        vector<double> values_in_time_model = params(time_model);

        Q_profile = vector<double>(pipe.profile.getPointCount(), values_in_time_model[3]); // задаем по трубе новый расход из временного ряда
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);

        quickest_ultimate_fv_solver solver_rho(*advection_model, density_wrapper);
        solver_rho.step(dt, values_in_time_model[0], values_in_time_model[0]);
        auto rho_profile = density_wrapper.current().vars;

        quickest_ultimate_fv_solver solver_nu(*advection_model, viscosity_wrapper);
        solver_nu.step(dt, values_in_time_model[1], values_in_time_model[1]);
        auto nu_profile = viscosity_wrapper.current().vars;
        
        pipe_model_PQ_cell_parties_t pipeModel(pipe, rho_profile, nu_profile, Q_profile[0], euler_direction);

        solve_euler<1>(pipeModel, euler_direction, values_in_time_model[2], &p_profile);

        std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
            [](double initial, double current) {return initial - current;  });

        buffer.advance(+1);
    } while (t < std::time(nullptr) + duration - dt);
}