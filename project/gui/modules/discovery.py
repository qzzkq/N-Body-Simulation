"""
Динамическое обнаружение алгоритмов и сценариев.

Алгоритмы — по содержимому папки `project/algo/`.
Сценарии  — парсинг `getName()` из `project/src/generators.cpp`.
"""

import os
import re


# Маппинг: stem файла в algo/ -> метаданные для GUI.
# Чтобы новый алгоритм появился в GUI:
#   1) положить <name>.cpp (или .cu) в project/algo/
#   2) добавить сюда строку с label и code (тот же код, что ждёт --algo в main.cpp)
ALGO_META = {
    "brut_force":      {"label": "Брутфорс  O(N²)",  "code": "1"},
    "barnes_hut":      {"label": "Барнс-Хат CPU",    "code": "2"},
    "barnes_hut_cuda": {"label": "Барнс-Хат CUDA",   "code": "3", "needs_cuda": True},
    "whfast":          {"label": "WHFast",            "code": "4"},
}


def discover_algorithms(project_root: str) -> list[dict]:
    """
    Возвращает список доступных алгоритмов в виде [{"label": str, "code": str}, ...]
    в порядке возрастания code. Если папки algo/ нет — пустой список (GUI покажет fallback).
    """
    algo_dir = os.path.join(project_root, "algo")
    if not os.path.isdir(algo_dir):
        return []

    found: list[dict] = []
    for fname in os.listdir(algo_dir):
        stem, ext = os.path.splitext(fname)
        if ext not in (".cpp", ".cu"):
            continue
        meta = ALGO_META.get(stem)
        if not meta:
            continue
        found.append({"label": meta["label"], "code": meta["code"]})

    found.sort(key=lambda d: int(d["code"]))
    return found


# Регексп ловит тело getName() — берём строковый литерал из return.
_GETNAME_RE = re.compile(
    r'std::string\s+getName\s*\(\s*\)\s*const\s+override\s*\{\s*return\s*"([^"]+)"\s*;\s*\}',
)


def discover_scenarios(project_root: str) -> list[str]:
    """
    Список имён сценариев в порядке их регистрации (он же индекс для --scenario).
    Источник — src/generators.cpp. Если файла нет или ничего не нашли — [].
    """
    src = os.path.join(project_root, "src", "generators.cpp")
    if not os.path.isfile(src):
        return []

    try:
        with open(src, "r", encoding="utf-8") as f:
            text = f.read()
    except OSError:
        return []

    return _GETNAME_RE.findall(text)
