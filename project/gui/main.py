"""
N-Body Simulation GUI
Дизайн: тёмный сайдбар + starfield-контент.
Запускать: ./launch.sh
"""

import sys
import os
import random

GUI_DIR = os.path.dirname(os.path.abspath(__file__))
if GUI_DIR not in sys.path:
    sys.path.insert(0, GUI_DIR)

import tkinter as tk
import customtkinter as ctk

from modules.config   import AppConfig
from modules.baker    import CalculatePage
from modules.replay   import PlayPage

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

PROJECT_ROOT = os.path.dirname(GUI_DIR)

# ── Цветовая палитра ────────────────────────────────────────
C_SIDEBAR      = "#1e2028"
C_SIDEBAR_BTN  = "#252830"
C_BTN_ACTIVE   = "#6b8fa5"
C_BG           = "#080b12"
C_TEXT         = "#ffffff"
C_TEXT_DIM     = "#8a8fa8"
C_BORDER       = "#2e3240"
C_START        = "#3d9de0"
C_START_HOVER  = "#5ab0f0"

SIDEBAR_W  = 270
FONT_TITLE = ("Georgia", 30, "bold")
FONT_SUB   = ("Georgia", 13)
FONT_BTN   = ("Helvetica Neue", 15, "bold")
FONT_HINT  = ("Helvetica Neue", 12)

PALETTE = dict(
    bg=C_BG, border=C_BORDER, text=C_TEXT,
    dim=C_TEXT_DIM, start=C_START, start_hover=C_START_HOVER,
    sidebar=C_SIDEBAR,
)


# ════════════════════════════════════════════════════════════
class StarfieldCanvas(tk.Canvas):
    """Звёздный фон."""
    def __init__(self, parent, **kw):
        super().__init__(parent, bg=C_BG, highlightthickness=0, **kw)
        self._stars: list = []
        self.bind("<Configure>", self._on_resize)

    def _on_resize(self, event):
        self.delete("all")
        self._stars.clear()
        n = int(event.width * event.height / 900)
        for _ in range(n):
            x = random.uniform(0, event.width)
            y = random.uniform(0, event.height)
            s = random.choice([0.6, 0.6, 0.8, 0.8, 1.0, 1.2, 1.8])
            a = random.uniform(0.3, 1.0)
            v = int(255 * a)
            col = f"#{v:02x}{v:02x}{v:02x}"
            self.create_oval(x - s, y - s, x + s, y + s,
                             fill=col, outline="", tags="star")


# ════════════════════════════════════════════════════════════
class SidebarButton(tk.Frame):
    """Кнопка навигации."""
    def __init__(self, parent, text: str, command=None, **kw):
        super().__init__(parent, bg=C_SIDEBAR_BTN,
                         highlightbackground=C_BORDER,
                         highlightthickness=1,
                         cursor="hand2", **kw)
        self._command = command
        self._active  = False
        self._lbl = tk.Label(self, text=text, bg=C_SIDEBAR_BTN,
                             fg=C_TEXT, font=FONT_BTN)
        self._lbl.pack(fill="both", expand=True, ipady=14)
        for w in (self, self._lbl):
            w.bind("<Button-1>", lambda _: self._command() if self._command else None)
            w.bind("<Enter>",    self._enter)
            w.bind("<Leave>",    self._leave)

    def _enter(self, _=None):
        if not self._active:
            c = "#2d3040"
            self.configure(bg=c); self._lbl.configure(bg=c)

    def _leave(self, _=None):
        c = C_BTN_ACTIVE if self._active else C_SIDEBAR_BTN
        self.configure(bg=c); self._lbl.configure(bg=c)

    def set_active(self, active: bool):
        self._active = active
        c = C_BTN_ACTIVE if active else C_SIDEBAR_BTN
        self.configure(bg=c); self._lbl.configure(bg=c)


# ════════════════════════════════════════════════════════════
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.config_obj = AppConfig(GUI_DIR)
        self.title("N-Body Simulation")
        self.geometry("1160x740")
        self.minsize(920, 600)
        self.configure(bg=C_SIDEBAR)
        self._current = None
        self._pages: dict[str, tk.Widget] = {}
        self._build()
        self._show("home")

    def _build(self):
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # ── Сайдбар ──────────────────────────────────────────
        sb = tk.Frame(self, bg=C_SIDEBAR, width=SIDEBAR_W)
        sb.grid(row=0, column=0, sticky="ns")
        sb.grid_propagate(False)

        tk.Label(sb, text="N-Body", bg=C_SIDEBAR,
                 fg=C_TEXT, font=FONT_TITLE).pack(
            anchor="w", padx=26, pady=(34, 4))
        tk.Label(sb, text="Your space lab", bg=C_SIDEBAR,
                 fg=C_TEXT_DIM, font=FONT_SUB).pack(
            anchor="w", padx=26, pady=(0, 44))

        self.btn_calc = SidebarButton(sb, "Calculate",
                                      command=lambda: self._show("calculate"))
        self.btn_calc.pack(fill="x")

        tk.Frame(sb, bg=C_BORDER, height=1).pack(fill="x")   # разделитель

        self.btn_play = SidebarButton(sb, "Play",
                                      command=lambda: self._show("play"))
        self.btn_play.pack(fill="x")

        # ── Контентная зона ───────────────────────────────────
        self._host = tk.Frame(self, bg=C_BG)
        self._host.grid(row=0, column=1, sticky="nsew")
        self._host.grid_rowconfigure(0, weight=1)
        self._host.grid_columnconfigure(0, weight=1)

        # Starfield — всегда позади
        self._stars = StarfieldCanvas(self._host)
        self._stars.place(x=0, y=0, relwidth=1, relheight=1)

        # ── Страницы ─────────────────────────────────────────
        self._pages["home"]      = self._make_home()
        self._pages["calculate"] = CalculatePage(self._host, PROJECT_ROOT,
                                                  self.config_obj, PALETTE)
        self._pages["play"]      = PlayPage(self._host, PROJECT_ROOT,
                                            self.config_obj, PALETTE)

    def _make_home(self) -> tk.Frame:
        f = tk.Frame(self._host, bg="")   # прозрачный — видны звёзды
        tk.Label(f, text="Choose the action", bg=C_BG,
                 fg=C_TEXT_DIM, font=FONT_HINT).place(
            relx=0.5, rely=0.97, anchor="s")
        return f

    def _show(self, page: str):
        if self._current:
            self._pages[self._current].place_forget()
        self._pages[page].place(x=0, y=0, relwidth=1, relheight=1)
        self._current = page
        self.btn_calc.set_active(page == "calculate")
        self.btn_play.set_active(page == "play")


if __name__ == "__main__":
    app = App()
    app.mainloop()
