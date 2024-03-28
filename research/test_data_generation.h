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

/// @brief Исходны данные и настроечные параметры
struct timeseries_generator_settings {
    std::time_t duration{ 300000 }; // Время моделирования, с
    std::time_t sample_time_min{200}; // Минимальное значение размаха шага, с
    std::time_t sample_time_max{400}; // Максимальное значение размаха шага, с
    double value_relative_decrement{ 0.0002 }; // Относительное минимальное отклонение значения параметров, доли
    double value_relative_increment{ 0.0002 }; // Относительное максимальное отклонение значения параметров, доли
    /// @brief Исходные данные, обязательно должны присутствовать два опциональных параметра
    /// @brief "rho_in" Плотность жидкости, (кг/м3)
    /// @brief "visc_in" Кинематическая вязкость, (м2/с)
    /// @brief "p_in" Давление на входе (опционально), (Па)
    /// @brief "p_out" Давление на выходе (опционально), (Па)
    /// @brief "Q" Расход по всей трубе (опционально), (м^3/с)
    vector<pair<string, double>> timeseries_initial_values = {
        { "rho_in", 860 },
        { "visc_in", 15e-6},
        { "p_in", 6e6},
        { "Q", 0.2 },
    };
};

/// @brief Класс для генерации синтетических временных рядов
class syntetic_time_series_generator {
public:
    /// @brief Конструктор класса
    /// @param settings Настройки генератора временных рядов
    syntetic_time_series_generator(const timeseries_generator_settings& settings)
        : settings_(settings),
        gen(rd()) {
        std::uniform_real_distribution<double> timeDis(settings_.sample_time_min, settings_.sample_time_max);

        for (const auto& param : settings_.timeseries_initial_values) {
            TimeVector timeValues;
            ParamVector paramValues;

            std::uniform_real_distribution<double> normalDis(param.second * (1 - settings_.value_relative_increment), param.second * (1 + settings_.value_relative_decrement));

            time_t timeStep = timeDis(gen);
            for (time_t time = std::time(nullptr); time <= std::time(nullptr) + settings_.duration; time += timeStep) {
                timeValues.push_back(time);
                timeStep = timeDis(gen);
            }

            std::transform(timeValues.begin(), timeValues.end(), std::back_inserter(paramValues),
                [&](time_t) { return normalDis(gen); });

            data.push_back({ timeValues, paramValues });
        }
    }
    /// @brief Применение скачка к временному ряду
    /// @param jump_time Время, когда происходит скачок
    /// @param jump_value Значение скачка
    /// @param paramName Имя параметра, к которому применяется скачок
    void apply_jump(time_t jump_time, double jump_value, const string& paramName) {
        for (size_t i = 0; i < settings_.timeseries_initial_values.size(); ++i) {
            if (settings_.timeseries_initial_values[i].first == paramName) {
                auto it = std::lower_bound(data[i].first.begin(), data[i].first.end(), time(nullptr) + jump_time);
                size_t position = std::distance(data[i].first.begin(), it);
                std::uniform_real_distribution<double> normalDis(jump_value * (1 - settings_.value_relative_increment), jump_value * (1 + settings_.value_relative_increment));
                for (size_t j = position; j < data[i].first.size(); ++j) {
                    double value = normalDis(gen);
                    data[i].second[j] = value;
                }
            }
        }
    }
    /// @brief Получение сгенерированных данных
    /// @return Вектор временных рядов
    vector<ParamPair> get_data() const {
        return data;
    }

private:
    timeseries_generator_settings settings_; // Настройки генератора
    vector<ParamPair> data;; // Данные временных рядов
    std::random_device rd; // Генератор случайных чисел
    std::mt19937 gen; // Генератор псевдослучайных чисел
};

TEST(Random, PrepareTimeSeries)
{
    timeseries_generator_settings settings;
    settings.duration = 350000; // Задаю время моделирования (опционально, по умолчанию 300000)
    settings.sample_time_min = 300; // Задаю минимальное значение размаха шага (опционально, по умолчанию 200)
    settings.sample_time_max = 450; // Задаю максимальное значение размаха шага (опционально, по умолчанию 400)
    settings.value_relative_decrement = 0.005; // Задаю относительное минимальное отклонение значения параметров (опционально, по умолчанию 0.0002)
    settings.value_relative_increment = 0.005; // Задаю относительное максимальное отклонение значения параметров (опционально, по умолчанию 0.0002)
    settings.timeseries_initial_values = {
        { "rho_in", 850 },
        { "visc_in", 17e-6},
        { "p_in", 5e6},
        { "Q", 0.3 },
    };
    syntetic_time_series_generator data_generator(settings);

    const time_t jump_time_rho = 100000;
    const double jump_value_rho = 870;
    data_generator.apply_jump(jump_time_rho, jump_value_rho, "rho_in");

    const time_t jump_time_Q = 150000;
    const double jump_value_Q = 0.1;
    data_generator.apply_jump(jump_time_Q, jump_value_Q, "Q");

    const auto data = data_generator.get_data();

    vector_timeseries_t params(data);

    // Задаём интересующий нас момент времени
    time_t test_time = static_cast<time_t>(std::time(nullptr) + 200000);

    // Интерополируем значения параметров в заданный момент времени
    vector<double> values_in_test_time = params(test_time);
}

#pragma once

/// @brief Тесты для солвера quickest_ultimate_fv_solver
class quick_with_quasi_stationary_model : public ::testing::Test {
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

TEST_F(quick_with_quasi_stationary_model, WorkingWithTimeSeries)
{
    // Объявляем структуру с исходными данными и настроечными параметрами
    timeseries_generator_settings settings;
    // Генерируем данные
    syntetic_time_series_generator data_time_series(settings); 
    // Получаем данные
    const auto data = data_time_series.get_data(); 
    // Помещаем временные ряды в вектор
    vector_timeseries_t params(data);

    const auto& x = advection_model->get_grid();
    double dx = x[1] - x[0]; // Шаг сетки

    ring_buffer_t<density_viscosity_cell_layer> buffer(2, pipe.profile.getPointCount());

    buffer[0].density = vector<double>(buffer[0].density.size(), 850); // Инициализация начальной плотности
    buffer[0].viscosity = vector<double>(buffer[0].viscosity.size(), 1e-5); // Инициализация начальной вязкости
    double p_initial = 6e6; // Давление на входе изначальное
    buffer.advance(+1);

    double v_max = 1; // Предполагаем скорость для Куранта = 1, скорость, больше чем во временных рядах и в профиле
    time_t dt = abs(dx/v_max); // Постоянный шаг по времени для Куранта = 1

    time_t t = std::time(nullptr); // Момент времени начала моделирования

    vector<double> initial_p_profile; // Изначальный профиль давлений
    vector<double> diff_p_profile = vector<double>(pipe.profile.getPointCount(), 0); // Дифференциальный профиль давлений

    do
    {
        vector<double> p_profile(pipe.profile.getPointCount()); // Профиль давлений

        auto density_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_density_quick_wrapper);

        auto viscosity_wrapper = buffer.get_buffer_wrapper(
            &density_viscosity_cell_layer::get_viscosity_quick_wrapper);

        int euler_direction = +1; // Задаем направление для Эйлера
        
        if (t == std::time(nullptr))
        {
            pipe_model_PQ_cell_parties_t pipeModel(pipe, density_wrapper.previous().vars, viscosity_wrapper.previous().vars, Q_profile[0], euler_direction);
            solve_euler<1>(pipeModel, euler_direction, p_initial, &p_profile);
            initial_p_profile = p_profile; // Получаем изначальный профиль давлений
        }
        t += dt;
        
        // Интерополируем значения параметров в заданный момент времени
        vector<double> values_in_time_model = params(t);

        Q_profile = vector<double>(pipe.profile.getPointCount(), values_in_time_model[3]); // задаем по трубе новый расход из временного ряда
        advection_model = std::make_unique<PipeQAdvection>(pipe, Q_profile);
        // Шаг по плотности
        quickest_ultimate_fv_solver solver_rho(*advection_model, density_wrapper);
        solver_rho.step(dt, values_in_time_model[0], values_in_time_model[0]);
        // Шаг по вязкости
        quickest_ultimate_fv_solver solver_nu(*advection_model, viscosity_wrapper);
        solver_nu.step(dt, values_in_time_model[1], values_in_time_model[1]);
        
        pipe_model_PQ_cell_parties_t pipeModel(pipe, density_wrapper.current().vars, viscosity_wrapper.current().vars, Q_profile[0], euler_direction);
        // Получаем новый профиль давлений
        solve_euler<1>(pipeModel, euler_direction, values_in_time_model[2], &p_profile);
        // Получаем дифференциальный профиль давлений
        std::transform(initial_p_profile.begin(), initial_p_profile.end(), p_profile.begin(), diff_p_profile.begin(),
            [](double initial, double current) {return initial - current;  });

        buffer.advance(+1);
    } while (t < std::time(nullptr) + settings.duration - dt);
}