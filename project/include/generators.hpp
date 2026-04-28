/**
 * @file generators.hpp
 * @brief Процедурная генерация начальных условий для N-body систем.
 *
 * Предоставляет расширяемую систему генераторов на основе паттерна
 * «Стратегия» (IScenario + ScenarioManager).
 *
 * Встроенные сценарии (регистрируются в CreateDefaultManager()):
 *  - **Кольцо из объектов** — тела равномерно распределены по кольцу
 *    на орбите вокруг центрального массивного тела.
 *  - **8 Кластеров** — тела разбиты на 8 кластеров, равномерно
 *    расположенных по орбите; каждый кластер орбитально устойчив.
 *
 * Новый сценарий добавляется за три шага:
 *  1. Создать класс, наследующий IScenario.
 *  2. Реализовать getName() и generate().
 *  3. Зарегистрировать через ScenarioManager::registerScenario() в
 *     CreateDefaultManager() внутри generators.cpp.
 */

#pragma once
#include <vector>
#include <string>
#include <memory>
#include <glm/glm.hpp>
#include "object.hpp"

/**
 * @struct GenParams
 * @brief Параметры генерации системы тел.
 *
 * Передаётся в IScenario::generate() и определяет масштаб и физические
 * характеристики генерируемой системы. Часть полей может игнорироваться
 * конкретным сценарием.
 */
struct GenParams {
    int    count;        ///< Количество тел (без учёта центрального тела, если оно есть).
    double centralMass;  ///< Масса центрального тела [M☉]. 0 — без центрального тела.
    double baseMass;     ///< Базовая масса рядовых тел [M☉].
    float  minRadius;    ///< Минимальный орбитальный радиус [AU].
    float  maxRadius;    ///< Максимальный орбитальный радиус [AU].
    float  spread;       ///< Разброс по Z-оси («толщина» диска) [AU].
    unsigned seed;       ///< Seed генератора случайных чисел. 0 — случайный.
};

/**
 * @class IScenario
 * @brief Абстрактный интерфейс сценария генерации начальных условий.
 *
 * Каждый наследник реализует конкретную конфигурацию системы тел.
 * Метод generate() должен:
 *  - Добавить тела в @p out (не очищать его перед добавлением).
 *  - Назначить корректные начальные скорости для орбитальной устойчивости
 *    (через формулу круговой орбиты v = √(G·M/r)).
 */
class IScenario {
public:
    virtual ~IScenario() = default;

    /**
     * @brief Возвращает человекочитаемое имя сценария.
     *
     * Имя отображается в меню выбора в main.cpp и в GUI.
     * @return Строка с именем (UTF-8, на русском).
     */
    virtual std::string getName() const = 0;

    /**
     * @brief Генерирует тела и добавляет их в массив.
     *
     * @param out    [out] Массив тел для добавления (append, не replace).
     * @param params Параметры генерации.
     */
    virtual void generate(std::vector<Object>& out, const GenParams& params) = 0;
};

/**
 * @class ScenarioManager
 * @brief Реестр сценариев генерации. Позволяет добавлять и запускать сценарии по индексу.
 *
 * Используется в main.cpp для отображения меню выбора сценария
 * и вызова нужного генератора.
 *
 * @par Пример:
 * @code
 * auto mgr = CreateDefaultManager();
 * for (auto& name : mgr->getNames())
 *     std::cout << name << "\n";
 * mgr->runScenario(0, objs, params);   // запустить первый сценарий
 * @endcode
 */
class ScenarioManager {
public:
    /**
     * @brief Регистрирует новый сценарий.
     *
     * Берёт владение объектом через unique_ptr.
     * @param scenario Указатель на реализацию IScenario.
     */
    void registerScenario(std::unique_ptr<IScenario> scenario);

    /**
     * @brief Возвращает список имён всех зарегистрированных сценариев.
     *
     * Порядок совпадает с порядком регистрации (и с индексами runScenario).
     * @return Вектор строк.
     */
    std::vector<std::string> getNames() const;

    /**
     * @brief Запускает сценарий по индексу.
     *
     * Вызывает IScenario::generate() с заданными параметрами.
     * При неверном индексе выводит ошибку в stderr.
     *
     * @param index  0-based индекс сценария.
     * @param out    [out] Целевой массив тел.
     * @param params Параметры генерации.
     */
    void runScenario(size_t index, std::vector<Object>& out, const GenParams& params);

    /**
     * @brief Проверяет валидность индекса сценария.
     * @param index Индекс для проверки.
     * @return @c true если индекс < количества зарегистрированных сценариев.
     */
    bool isValidIndex(size_t index) const;

private:
    std::vector<std::unique_ptr<IScenario>> scenarios_; ///< Список зарегистрированных сценариев.
};

/**
 * @brief Создаёт ScenarioManager с предустановленными сценариями.
 *
 * Регистрирует:
 *  - RingScenario    («Кольцо из объектов»)
 *  - ClusterScenario («8 Кластеров»)
 *
 * Вызывается один раз в main() для инициализации меню выбора.
 *
 * @return unique_ptr на настроенный ScenarioManager.
 */
std::unique_ptr<ScenarioManager> CreateDefaultManager();