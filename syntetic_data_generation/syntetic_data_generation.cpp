#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>

using namespace std;

// Типы синонимов для улучшения читаемости кода
using TimeVector = vector<time_t>;
using ParamVector = vector<double>;
using ParamPair = pair<TimeVector, ParamVector>;

class SyntheticDataGenerator {
public:
    SyntheticDataGenerator(const vector<pair<string, double>>& params)
        : parameters(params) {}

    void createDataRows(double duration) {
        data.clear();

        // Генерация случайных чисел

        uniform_real_distribution<double> timeDis(200, 400); // Для времени

        // Генерация данных для каждого параметра
        for (const auto& param : parameters) {
            TimeVector timeValues;
            ParamVector paramValues;

            uniform_real_distribution<double> normalDis(param.second * 0.9998, param.second * 1.0002); // Для значений

            // Генерация времени с постоянным шагом
            double timeStep = timeDis(gen);
            for (double time = std::time(nullptr); time <= std::time(nullptr) + duration; time += timeStep) {
                timeValues.push_back(static_cast<time_t>(time));
                timeStep = timeDis(gen); // Переменный шаг
            }

            // Генерация значений с разбросом
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
                // Генерация значений с разбросом
                uniform_real_distribution<double> normalDis(jumpValue * 0.9998, jumpValue * 1.0002); // Для значений
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
};

int main() {
    const vector<pair<string, double>> parameters = {
        { "rho_in", 860 },
        { "visc_in", 15e-6},
        { "p_in", 6e6},
        { "Q", 0.2 },
    };

    const double duration = 350000;

    SyntheticDataGenerator dataGenerator(parameters);
    dataGenerator.createDataRows(duration);

    const double jumpTime = 100000;
    const double jumpValue = 870;
    dataGenerator.applyJump(jumpTime, jumpValue, "rho_in");

    const auto data = dataGenerator.getData();

    return 0;
}