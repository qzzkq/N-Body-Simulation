"""
Play Page — starfield-фон, по центру «Put the path» + поле ввода,
кнопка Start внизу по центру.
Точно по макету из PDF (страница 2).
"""

import os
import subprocess
from tkinter import filedialog
import tkinter as tk
import customtkinter as ctk

FONT_HEADER  = ("Helvetica Neue", 15, "bold")
FONT_LABEL   = ("Helvetica Neue", 11)
FONT_START   = ("Helvetica Neue", 15, "bold")
FONT_META    = ("Courier New",    10)
FONT_SMALL   = ("Helvetica Neue", 10)


class RoundedButton(tk.Canvas):
    def __init__(self, parent, text, command=None, width=200, height=46,
                 bg_color="#3d9de0", bg_hover="#5ab0f0",
                 text_color="#ffffff", font=None, **kw):
        if "bg" not in kw:
            kw["bg"] = parent["bg"]
        super().__init__(parent, width=width, height=height,
                         highlightthickness=0, bd=0, relief="flat", cursor="hand2", **kw)
        self._text     = text
        self._command  = command
        self._bg       = bg_color
        self._bg_hover = bg_hover
        self._tc       = text_color
        self._font     = font or FONT_START
        self._hovered  = False
        self.bind("<Configure>", self._draw)
        self.bind("<Enter>",     lambda _: self._hover(True))
        self.bind("<Leave>",     lambda _: self._hover(False))
        self.bind("<Button-1>",  lambda _: command() if command else None)
        self._draw()

    def _draw(self, _=None):
        self.delete("all")
        w = self.winfo_width()  or int(self["width"])
        h = self.winfo_height() or int(self["height"])
        r   = h // 2
        col = self._bg_hover if self._hovered else self._bg
        self.create_arc(0, 0, r*2, h,     start=90,  extent=180, fill=col, outline="")
        self.create_arc(w-r*2, 0, w, h,   start=270, extent=180, fill=col, outline="")
        self.create_rectangle(r, 0, w-r, h, fill=col, outline="")
        self.create_text(w//2, h//2, text=self._text, fill=self._tc, font=self._font, anchor="center")

    def _hover(self, state: bool):
        self._hovered = state
        self._draw()


class RoundedEntry(tk.Canvas):
    """Поле ввода со скруглёнными краями — как в макете."""
    def __init__(self, parent, palette: dict, width=520, height=52, **kw):
        if "bg" not in kw:
            kw["bg"] = parent["bg"]
        super().__init__(parent, width=width, height=height,
                         highlightthickness=0, bd=0, relief="flat", **kw)
        self._pal = palette
        self._bg_c = "#12141e"
        self.bind("<Configure>", self._draw_bg)
        self._var = tk.StringVar()
        self._inner = tk.Entry(
            self, textvariable=self._var,
            bg=self._bg_c, fg=palette["text"],
            insertbackground=palette["text"],
            relief="flat", font=("Helvetica Neue", 12),
            justify="center",
        )
        self._win_id = None
        self.bind("<Configure>", self._on_cfg)
        self._draw_bg()

    def _on_cfg(self, event):
        self._draw_bg()
        w, h = event.width, event.height
        r = h // 2
        if self._win_id:
            self.delete(self._win_id)
        self._win_id = self.create_window(
            r + 4, h // 2,
            width=w - r * 2 - 8, height=h - 16,
            window=self._inner, anchor="w",
        )

    def _draw_bg(self, _=None):
        self.delete("bg")
        w = self.winfo_width()  or int(self["width"])
        h = self.winfo_height() or int(self["height"])
        r = h // 2
        col = self._bg_c
        bdr = self._pal["border"]
        self.create_arc(0, 0, r*2, h,   start=90,  extent=180, fill=bdr, outline="", tags="bg")
        self.create_arc(w-r*2, 0, w, h, start=270, extent=180, fill=bdr, outline="", tags="bg")
        self.create_rectangle(r, 0, w-r, h, fill=bdr, outline="", tags="bg")
        m = 1
        self.create_arc(m, m, r*2-m, h-m,       start=90,  extent=180, fill=col, outline="", tags="bg")
        self.create_arc(w-r*2+m, m, w-m, h-m,   start=270, extent=180, fill=col, outline="", tags="bg")
        self.create_rectangle(r, m, w-r, h-m,   fill=col, outline="", tags="bg")

    def get(self) -> str:
        return self._var.get()

    def set(self, v: str):
        self._var.set(v)

    def bind_entry(self, event, func):
        self._inner.bind(event, func)


class PlayPage(tk.Frame):
    def __init__(self, parent, project_root: str, config, palette: dict):
        super().__init__(parent, bg=palette["bg"])
        self.project_root = project_root
        self.config       = config
        self.pal          = palette
        self._meta: dict  = {}
        self._build()

    # ════════════════════════════════════════════════════════
    def _build(self):
        center = tk.Frame(self, bg=self.pal["bg"])
        center.place(relx=0.5, rely=0.42, anchor="center")

        # ── «Put the path» ────────────────────────────────
        ctk.CTkLabel(
            center, 
            text="Put the path",
            fg_color=self.pal["border"],   
            text_color=self.pal["text"],   
            font=ctk.CTkFont(family="Helvetica Neue", size=14, weight="bold"),
            width=300, 
            height=46, 
            corner_radius=23              
        ).pack(pady=(0, 6))

        # ── Поле ввода ────────────────────────────────────
        entry_row = tk.Frame(center, bg=self.pal["bg"])
        entry_row.pack(pady=(6, 0))

        self.path_entry = RoundedEntry(entry_row, self.pal,
                                       width=500, height=50,
                                       bg=self.pal["bg"])
        self.path_entry.pack(side="left")

        browse_btn = tk.Button(
            entry_row, text="📂",
            bg=self.pal["border"], fg=self.pal["text"],
            relief="flat", cursor="hand2", padx=10,
            activebackground="#3a3f55",
            font=("Helvetica Neue", 16),
            command=self._browse,
        )
        browse_btn.pack(side="left", padx=(8, 0), ipady=8)

        # ── Метаданные (маленький блок под полем) ─────────
        self.meta_frame = tk.Frame(center, bg=self.pal["bg"])
        self.meta_frame.pack(pady=(10, 0))
        self.meta_lbl = tk.Label(
            self.meta_frame, text="",
            bg=self.pal["bg"], fg=self.pal["dim"],
            font=FONT_META, justify="center",
        )
        self.meta_lbl.pack()

        ctk.CTkButton(
            center, text="🔍 Read metadata", height=28, width=160,
            fg_color=self.pal["border"], hover_color="#3a3f55",
            font=ctk.CTkFont(size=11),
            command=self._read_meta,
        ).pack(pady=(8, 0))

        # ── Start — внизу по центру ───────────────────────
        self.start_btn = RoundedButton(
            self, text="Start", width=220, height=48,
            bg_color=self.pal["start"],
            bg_hover=self.pal["start_hover"],
            command=self._on_start,
        )
        self.start_btn.place(relx=0.5, rely=0.93, anchor="center")

    # ════════════════════════════════════════════════════════
    def _draw_pill(self, canvas: tk.Canvas, text: str,
                   bg: str, text_color: str, font):
        canvas.update_idletasks()
        w = canvas.winfo_width()  or int(canvas["width"])
        h = canvas.winfo_height() or int(canvas["height"])
        r = h // 2
        canvas.create_arc(0, 0, r*2, h, start=90, extent=180,
                          fill=bg, outline="")
        canvas.create_arc(w-r*2, 0, w, h, start=270, extent=180,
                          fill=bg, outline="")
        canvas.create_rectangle(r, 0, w-r, h, fill=bg, outline="")
        canvas.create_text(w//2, h//2, text=text,
                           fill=text_color, font=font, anchor="center")

    # ── Browse ────────────────────────────────────────────
    def _browse(self):
        data = os.path.join(self.project_root,
                            self.config.get("data_dir", "data"))
        path = filedialog.askopenfilename(
            initialdir=data if os.path.isdir(data) else ".",
            filetypes=[("HDF5 files", "*.h5 *.hdf5"), ("All files", "*.*")],
        )
        if path:
            self.path_entry.set(path)

    # ── Чтение метаданных ─────────────────────────────────
    def _read_meta(self):
        path = self.path_entry.get().strip()
        if not path:
            self.meta_lbl.configure(text="No path specified")
            return
        full = path if os.path.isabs(path) \
               else os.path.join(self.project_root, path)
        if not os.path.isfile(full):
            self.meta_lbl.configure(text=f"File not found: {path}")
            return
        try:
            import h5py
            sz = os.path.getsize(full) / 1024 / 1024
            with h5py.File(full, "r") as f:
                bodies = int(f.attrs.get("num_bodies", 0))
                frames = int(f.attrs.get("num_frames", 0))
                dt     = float(f.attrs.get("dt", 0.0))
            dur = frames * dt
            self._meta = dict(bodies=bodies, frames=frames, dt=dt, duration=dur)
            self.meta_lbl.configure(
                text=f"Bodies: {bodies}   Frames: {frames}   "
                     f"dt: {dt:.6f}   Duration: {dur:.2f} yr   "
                     f"Size: {sz:.2f} MB")
        except ImportError:
            self.meta_lbl.configure(text="h5py not installed")
        except Exception as e:
            self.meta_lbl.configure(text=f"Error: {e}")

    # ── Start ─────────────────────────────────────────────
    def _on_start(self):
        path = self.path_entry.get().strip().replace('"', '') 
        
        if not path:
            self.meta_lbl.configure(text="Сначала выберите .h5 файл")
            return
            
        if not os.path.exists(path):
            self.meta_lbl.configure(text="Ошибка: Файл не найден!", fg="#f87171")
            return
            
        binary_name = self.config.get("replay_bin", "replay")
        binary_path = os.path.join(self.project_root, binary_name)
        
        if not os.path.exists(binary_path):
            self.meta_lbl.configure(text=f"Ошибка: Исполняемый файл '{binary_name}' не найден.", fg="#f87171")
            return
            
        self.meta_lbl.configure(text="Запуск Replay...")
        
        # replay.cpp ожидает путь к файлу первым аргументом без флагов
        cmd = [binary_path, path]
        subprocess.Popen(cmd, cwd=self.project_root)
