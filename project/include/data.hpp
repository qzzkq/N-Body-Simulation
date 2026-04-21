/**
 * @file data.hpp
 * @brief Ввод/вывод данных симуляции: HDF5-файлы и текстовые конфигурации.
 *
 * Модуль реализует два формата хранения:
 *
 * **1. Формат начальных условий (`.h5` начальные данные / `.txt`)**
 *  Используется для задания начального состояния системы.
 *  В HDF5 хранится датасет @c Particles со структурой Particle.
 *  В TXT — текстовый формат: по одной строке на тело.
 *
 * **2. Формат записи симуляции (`frames.h5`)**
 *  Создаётся в процессе работы. Содержит:
 *  - Атрибуты файла: @c num_bodies, @c num_frames, @c dt.
 *  - Датасет @c initial_data — начальные позиции и массы.
 *  - Датасет @c tracks — 2D-массив [кадр × тело] позиций (float).
 *  - Датасет @c body_colors — цвета тел (float RGBA).
 *
 *  Replay читает именно этот формат.
 */

#ifndef DATA_HPP
#define DATA_HPP
#pragma once

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "H5Cpp.h"
#include "object.hpp"
#include "graphic_state.hpp"

/**
 * @struct Particle
 * @brief Запись одного тела в HDF5-датасете начальных условий.
 *
 * Используется функциями Writer() / Reader() для сериализации
 * начального состояния системы в HDF5 составном типе (CompType).
 * Радиус не хранится — восстанавливается из mass и density.
 */
struct Particle {
    glm::dvec3 position;                      ///< Позиция [AU].
    glm::dvec3 velocity;                      ///< Скорость [AU/год].
    double     mass;                          ///< Масса [M☉].
    double     radius;                        ///< Радиус [AU].
    glm::vec4  color{1.0f, 1.0f, 1.0f, 1.0f};///< RGBA-цвет (опциональный, float).
};

// ── API записи симуляции (frames.h5) ────────────────────────────────────────

/**
 * @brief Создаёт новый HDF5-файл для записи траекторий.
 *
 * Инициализирует структуру файла:
 *  - Атрибуты: @c num_bodies, @c dt.
 *  - Датасет @c initial_data (позиции, массы, радиусы тел).
 *  - Датасет @c tracks (расширяемый: [кадр × тело], float позиции).
 *  - Датасет @c body_colors (RGBA цвета тел).
 *
 * @param fileName   Путь к создаваемому файлу (например, @c "data/frames.h5").
 * @param numBodies  Количество тел в системе.
 * @param dt         Симуляционный интервал между сохраняемыми кадрами
 *                   = fixedDt × saveIntervalSteps [год].
 * @return Открытый HDF5-файл. Вызывающая сторона отвечает за закрытие через CloseSimulationFile().
 */
H5::H5File CreateSimulationFile(const std::string& fileName, std::size_t numBodies, double dt);

/**
 * @brief Записывает один кадр (снимок позиций всех тел) в файл.
 *
 * Дозаписывает строку в расширяемый датасет @c tracks.
 * Позиции хранятся как @c float (не double) для экономии места.
 *
 * @param file        Открытый HDF5-файл (из CreateSimulationFile()).
 * @param objs        Текущее состояние тел.
 * @param graphics    Визуальные состояния (для записи цветов при первом кадре).
 * @param frameIndex  Индекс текущего кадра (0-based, используется как offset в датасете).
 */
void WriteSimulationFrame(H5::H5File& file,
                          const std::vector<Object>& objs,
                          const std::vector<GraphicState>& graphics,
                          std::size_t frameIndex);

/**
 * @brief Завершает запись файла симуляции.
 *
 * Обновляет атрибут @c num_frames и закрывает файл HDF5.
 *
 * @param file         Открытый HDF5-файл.
 * @param totalFrames  Финальное число записанных кадров.
 */
void CloseSimulationFile(H5::H5File& file, std::size_t totalFrames);

// ── API начальных условий (legacy) ──────────────────────────────────────────

/**
 * @brief Записывает начальные условия системы в HDF5-датасет.
 *
 * Сохраняет вектор Particle в датасет @p dsetName файла @p fileName.
 * Используется DataCreator для создания начальных .h5 файлов.
 *
 * @param fileName   Путь к HDF5-файлу.
 * @param dsetName   Имя датасета (например, @c "Particles").
 * @param particles  Данные для записи.
 */
void Writer(const std::string& fileName,
            const std::string& dsetName,
            const std::vector<Particle>& particles);

/**
 * @brief Читает начальные условия из HDF5-датасета.
 *
 * @param fileName             Путь к HDF5-файлу.
 * @param dsetName             Имя датасета.
 * @param outFileHadColorMember  Если не @c nullptr — устанавливается в @c true,
 *                             если датасет содержит поле @c color.
 * @return Вектор Particle с загруженными данными.
 */
std::vector<Particle> Reader(const std::string& fileName,
                             const std::string& dsetName,
                             bool* outFileHadColorMember = nullptr);

// ── Высокоуровневые функции загрузки ────────────────────────────────────────

/**
 * @brief Возвращает список .h5 файлов в директории.
 *
 * @param dir Директория для поиска (по умолчанию @c "data").
 * @return Отсортированный вектор путей (например, @c "data/solar_system.h5").
 */
std::vector<std::string> ListH5Files(const std::string& dir = "data");

/**
 * @brief Загружает начальные условия из .h5 файла в массивы Object и GraphicState.
 *
 * Вызывает Reader() и конвертирует Particle → Object.
 * Если файл содержит поле color — заполняет @p outGraphics.
 * После загрузки вызывает physics::colorFromMass() если цвета отсутствуют.
 *
 * @param filePath          Путь к HDF5-файлу с начальными условиями.
 * @param dsetName          Имя датасета (обычно @c "Particles").
 * @param outObjs           [out] Загруженные тела.
 * @param outGraphics       [out] Визуальные состояния (может быть @c nullptr).
 * @param outHadColorInFile [out] @c true если файл содержал цвета (может быть @c nullptr).
 * @return @c true при успешной загрузке.
 */
bool LoadObjectsFromFile(const std::string& filePath,
                         const std::string& dsetName,
                         std::vector<Object>& outObjs,
                         std::vector<GraphicState>* outGraphics = nullptr,
                         bool* outHadColorInFile = nullptr);

/**
 * @brief Загружает начальные условия из текстового файла.
 *
 * Формат файла: одна строка = одно тело.
 * Каждая строка: @c name x y z vx vy vz mass density r g b
 * Поля разделяются пробелами; строки начинающиеся с @c # — комментарии.
 *
 * @param filePath    Путь к .txt файлу.
 * @param outObjs     [out] Загруженные тела.
 * @param outGraphics [out] Визуальные состояния с цветами из файла (может быть @c nullptr).
 * @return @c true при успешной загрузке хотя бы одного тела.
 */
bool LoadSystemFromTextFile(const std::string& filePath,
                            std::vector<Object>& outObjs,
                            std::vector<GraphicState>* outGraphics = nullptr);

/**
 * @brief Сохраняет финальное состояние системы в текстовый файл.
 *
 * Используется для сохранения конечного состояния симуляции при выходе
 * (в том числе по Ctrl+C). Файл совместим с LoadSystemFromTextFile().
 *
 * @param filePath  Путь к сохраняемому файлу (например, @c "data/frames_final.txt").
 * @param objs      Тела для сохранения.
 * @return @c true при успехе.
 */
bool SaveSystemToTextFile(const std::string& filePath,
                          const std::vector<Object>& objs);

#endif // DATA_HPP

