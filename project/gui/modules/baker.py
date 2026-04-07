"""
Calculate Page — 4 панели на starfield-фоне + прогресс-бар + кнопка Start.

Панели (2×2):
  [0] Алгоритм + Render mode
  [1] Источник начальных условий
  [2] Физика (dt, save interval)
  [3] Вывод (filename, режим)
"""

import os
from tkinter import filedialog
import tkinter as tk
import customtkinter as ctk

# ── Данные из main.cpp ─────────────────────────────────────
RENDER_MODES   = [("Сферы", "1"), ("Кубы", "0"), ("Точки", "2")]
ALGORITHMS     = [("Брутфорс  O(N²)", "1"), ("Барнс-Хат CPU", "2"), ("Барнс-Хат CUDA", "3")]
SOURCES        = ["HDF5 (.h5)", "TXT (.txt)", "Случайная генерация"]
SCENARIOS      = ["Кольцо из объектов", "8 Кластеров"]

FONT_PANEL_HDR = ("Helvetica Neue", 13, "bold")
FONT_LABEL     = ("Helvetica Neue", 11)
FONT_SMALL     = ("Helvetica Neue", 10)
FONT_START     = ("Helvetica Neue", 15, "bold")
FONT_STATUS    = ("Helvetica Neue", 11)
FONT_PROGRESS  = ("Helvetica Neue", 10)


class RoundedButton(tk.Canvas):
    """Таблетообразная кнопка как в макете."""
    def __init__(self, parent, text, command=None, width=130, height=42,
                 bg_color="#3d9de0", bg_hover="#5ab0f0",
                 text_color="#ffffff", font=None, **kw):
        super().__init__(parent, width=width, height=height,
                         highlightthickness=0, cursor="hand2", **kw)
        self._text      = text
        self._command   = command
        self._bg        = bg_color
        self._bg_hover  = bg_hover
        self._tc        = text_color
        self._font      = font or FONT_START
        self._hovered   = False
        self.bind("<Configure>", self._draw)
        self.bind("<Enter>",     lambda _: self._hover(True))
        self.bind("<Leave>",     lambda _: self._hover(False))
        self.bind("<Button-1>",  lambda _: command() if command else None)
        self._draw()

    def _draw(self, _=None):
        self.delete("all")
        w, h = self.winfo_width() or int(self["width"]), \
               self.winfo_height() or int(self["height"])
        r   = h // 2
        col = self._bg_hover if self._hovered else self._bg
        # рисуем скруглённый прямоугольник
        self.create_arc(0, 0, r*2, h, start=90, extent=180, fill=col, outline="")
        self.create_arc(w - r*2, 0, w, h, start=270, extent=180, fill=col, outline="")
        self.create_rectangle(r, 0, w - r, h, fill=col, outline="")
        self.create_text(w//2, h//2, text=self._text,
                         fill=self._tc, font=self._font)

    def _hover(self, state: bool):
        self._hovered = state
        self._draw()


class Panel(tk.Frame):
    """Тёмная панель с тонкой рамкой — как в макете."""
    def __init__(self, parent, palette: dict, **kw):
        super().__init__(parent,
                         bg=palette["bg"],
                         highlightbackground=palette["border"],
                         highlightthickness=1,
                         **kw)
        self._pal = palette

    def section(self, text: str) -> tk.Label:
        lbl = tk.Label(self, text=text, bg=self._pal["bg"],
                       fg=self._pal["text"], font=FONT_PANEL_HDR, anchor="w")
        lbl.pack(fill="x", padx=14, pady=(14, 6))
        return lbl

    def label(self, text: str) -> tk.Label:
        lbl = tk.Label(self, text=text, bg=self._pal["bg"],
                       fg=self._pal["dim"], font=FONT_LABEL, anchor="w")
        lbl.pack(fill="x", padx=14, pady=(4, 0))
        return lbl

    def radio_group(self, var: tk.Variable, options: list[tuple[str, str]]):
        for label, val in options:
            rb = tk.Radiobutton(
                self, text=label, variable=var, value=val,
                bg=self._pal["bg"], fg=self._pal["text"],
                selectcolor=self._pal["bg"],
                activebackground=self._pal["bg"],
                activeforeground=self._pal["text"],
                font=FONT_LABEL, anchor="w",
            )
            rb.pack(fill="x", padx=24, pady=1)

    def entry(self, default: str = "", placeholder: str = "") -> tk.Entry:
        e = tk.Entry(self,
                     bg="#12141e", fg=self._pal["text"],
                     insertbackground=self._pal["text"],
                     relief="flat",
                     highlightbackground=self._pal["border"],
                     highlightthickness=1,
                     font=FONT_LABEL)
        if default:
            e.insert(0, default)
        elif placeholder:
            e.insert(0, placeholder)
            e.config(fg=self._pal["dim"])
            def _focus_in(ev, en=e, ph=placeholder):
                if en.get() == ph:
                    en.delete(0, "end")
                    en.config(fg=self._pal["text"])
            def _focus_out(ev, en=e, ph=placeholder):
                if not en.get():
                    en.insert(0, ph)
                    en.config(fg=self._pal["dim"])
            e.bind("<FocusIn>",  _focus_in)
            e.bind("<FocusOut>", _focus_out)
        e.pack(fill="x", padx=14, pady=(2, 4), ipady=4)
        return e

    def separator(self):
        tk.Frame(self, bg=self._pal["border"], height=1).pack(
            fill="x", padx=14, pady=6)

    def combo(self, values: list[str], var: tk.StringVar) -> ctk.CTkComboBox:
        cb = ctk.CTkComboBox(self, values=values, variable=var,
                             fg_color="#12141e", border_color=self._pal["border"],
                             button_color=self._pal["border"],
                             text_color=self._pal["text"],
                             font=ctk.CTkFont(size=11),
                             dropdown_fg_color="#1a1c26")
        cb.pack(fill="x", padx=14, pady=(2, 4))
        return cb


class CalculatePage(tk.Frame):
    def __init__(self, parent, project_root: str, config, palette: dict):
        super().__init__(parent, bg=palette["bg"])
        self.project_root = project_root
        self.config       = config
        self.pal          = palette
        self.data_dir     = os.path.join(project_root, config.get("data_dir", "data"))
        self._progress    = 0.0

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self._build()

    # ════════════════════════════════════════════════════════
    def _build(self):
        # ── Сетка 2×2 ────────────────────────────────────────
        grid = tk.Frame(self, bg=self.pal["bg"])
        grid.place(relx=0.03, rely=0.04, relwidth=0.94, relheight=0.86)
        grid.grid_rowconfigure(0, weight=1)
        grid.grid_rowconfigure(1, weight=1)
        grid.grid_columnconfigure(0, weight=1)
        grid.grid_columnconfigure(1, weight=1)

        GAP = 10
        p0 = Panel(grid, self.pal)
        p1 = Panel(grid, self.pal)
        p2 = Panel(grid, self.pal)
        p3 = Panel(grid, self.pal)

        p0.grid(row=0, column=0, sticky="nsew", padx=(0, GAP//2), pady=(0, GAP//2))
        p1.grid(row=0, column=1, sticky="nsew", padx=(GAP//2, 0), pady=(0, GAP//2))
        p2.grid(row=1, column=0, sticky="nsew", padx=(0, GAP//2), pady=(GAP//2, 0))
        p3.grid(row=1, column=1, sticky="nsew", padx=(GAP//2, 0), pady=(GAP//2, 0))

        self._build_panel_algo(p0)
        self._build_panel_source(p1)
        self._build_panel_physics(p2)
        self._build_panel_output(p3)

        # ── Нижняя полоса ─────────────────────────────────────
        self._build_bottom_bar()

    # ──────────────────────────────────────────────────────
    def _build_panel_algo(self, p: Panel):
        p.section("Algorithm")
        self.algo_var = tk.StringVar(value="2")
        p.radio_group(self.algo_var, ALGORITHMS)
        p.separator()
        p.section("Render mode")
        self.render_var = tk.StringVar(value="1")
        p.radio_group(self.render_var, RENDER_MODES)

    def _build_panel_source(self, p: Panel):
        p.section("Initial conditions")

        self.source_var = tk.StringVar(value=SOURCES[0])

        seg_frame = tk.Frame(p, bg=self.pal["bg"])
        seg_frame.pack(fill="x", padx=14, pady=(0, 8))
        self._src_btns: list[RoundedButton] = []
        for i, s in enumerate(SOURCES):
            lbl = s.split("(")[0].strip()
            btn = RoundedButton(
                seg_frame, text=lbl, width=88, height=28,
                bg_color=self.pal["border"],
                bg_hover="#3a3f55",
                text_color=self.pal["text"],
                font=("Helvetica Neue", 10),
                command=lambda v=s, i=i: self._set_source(v, i),
            )
            btn.grid(row=0, column=i, padx=2)
        seg_frame.grid_columnconfigure(0, weight=1)
        seg_frame.grid_columnconfigure(1, weight=1)
        seg_frame.grid_columnconfigure(2, weight=1)

        # sub-frames
        self._sf_hdf5 = tk.Frame(p, bg=self.pal["bg"])
        self._sf_txt  = tk.Frame(p, bg=self.pal["bg"])
        self._sf_rand = tk.Frame(p, bg=self.pal["bg"])

        self._build_sf_hdf5(self._sf_hdf5)
        self._build_sf_txt(self._sf_txt)
        self._build_sf_rand(self._sf_rand)

        self._set_source(SOURCES[0], 0)

    def _build_sf_hdf5(self, f):
        tk.Label(f, text="Select .h5 file", bg=self.pal["bg"],
                 fg=self.pal["dim"], font=FONT_LABEL).pack(
            fill="x", padx=14, pady=(4, 0))
        h5list = self._list_h5()
        self.h5_var   = tk.StringVar(value=h5list[0])
        cb = ctk.CTkComboBox(f, values=h5list, variable=self.h5_var,
                             fg_color="#12141e",
                             border_color=self.pal["border"],
                             button_color=self.pal["border"],
                             text_color=self.pal["text"],
                             font=ctk.CTkFont(size=11),
                             dropdown_fg_color="#1a1c26")
        cb.pack(fill="x", padx=14, pady=(4, 2))
        ctk.CTkButton(f, text="↻ Refresh", height=24, width=100,
                      fg_color=self.pal["border"],
                      hover_color="#3a3f55",
                      font=ctk.CTkFont(size=11),
                      command=lambda: cb.configure(values=self._list_h5())
                      ).pack(anchor="w", padx=14, pady=2)

    def _build_sf_txt(self, f):
        tk.Label(f, text="TXT file path", bg=self.pal["bg"],
                 fg=self.pal["dim"], font=FONT_LABEL).pack(
            fill="x", padx=14, pady=(4, 0))
        row = tk.Frame(f, bg=self.pal["bg"])
        row.pack(fill="x", padx=14, pady=(4, 2))
        row.grid_columnconfigure(0, weight=1)
        self.txt_entry = tk.Entry(row, bg="#12141e", fg=self.pal["text"],
                                  insertbackground=self.pal["text"],
                                  relief="flat",
                                  highlightbackground=self.pal["border"],
                                  highlightthickness=1,
                                  font=FONT_LABEL)
        self.txt_entry.grid(row=0, column=0, sticky="ew", ipady=4)
        tk.Button(row, text="📂", bg=self.pal["border"], fg=self.pal["text"],
                  relief="flat", cursor="hand2", padx=6,
                  command=self._browse_txt).grid(row=0, column=1, padx=(4, 0))

    def _build_sf_rand(self, f):
        tk.Label(f, text="Body count", bg=self.pal["bg"],
                 fg=self.pal["dim"], font=FONT_LABEL).pack(
            fill="x", padx=14, pady=(4, 0))
        self.rand_count = tk.Entry(f, bg="#12141e", fg=self.pal["text"],
                                   insertbackground=self.pal["text"],
                                   relief="flat",
                                   highlightbackground=self.pal["border"],
                                   highlightthickness=1,
                                   font=FONT_LABEL)
        self.rand_count.insert(0, self.config.get("default_bodies", "100"))
        self.rand_count.pack(fill="x", padx=14, pady=(4, 4), ipady=4)

        tk.Label(f, text="Scenario", bg=self.pal["bg"],
                 fg=self.pal["dim"], font=FONT_LABEL).pack(
            fill="x", padx=14, pady=(4, 0))
        self.scenario_var = tk.StringVar(value="1")
        for i, name in enumerate(SCENARIOS):
            tk.Radiobutton(f, text=name, variable=self.scenario_var,
                           value=str(i + 1),
                           bg=self.pal["bg"], fg=self.pal["text"],
                           selectcolor=self.pal["bg"],
                           activebackground=self.pal["bg"],
                           font=FONT_LABEL).pack(
                fill="x", padx=24, pady=1)

    def _set_source(self, value: str, idx: int):
        self.source_var.set(value)
        for f in (self._sf_hdf5, self._sf_txt, self._sf_rand):
            f.pack_forget()
        [self._sf_hdf5, self._sf_txt, self._sf_rand][idx].pack(fill="x")

    def _build_panel_physics(self, p: Panel):
        p.section("Physics")
        p.label("Time step  dt  (years/step)")
        self.dt_entry = p.entry(
            default=self.config.get("default_dt", "0.000273785"))
        p.separator()
        p.label("Save every N steps")
        self.save_n_entry = p.entry(
            default=self.config.get("default_save_n", "10"))
        p.separator()
        p.section("Mode")
        self.rt_var = tk.BooleanVar(value=True)

        rt_row = tk.Frame(p, bg=self.pal["bg"])
        rt_row.pack(fill="x", padx=14, pady=(4, 0))
        tk.Label(rt_row, text="Real-time", bg=self.pal["bg"],
                 fg=self.pal["text"], font=FONT_LABEL).grid(
            row=0, column=0, sticky="w")
        self._rt_switch = ctk.CTkSwitch(
            rt_row, text="", variable=self.rt_var,
            command=self._on_rt_toggle,
            onvalue=True, offvalue=False,
            progress_color=self.pal["start"],
        )
        self._rt_switch.grid(row=0, column=1, padx=12)
        tk.Label(rt_row, text="Bake", bg=self.pal["bg"],
                 fg=self.pal["dim"], font=FONT_LABEL).grid(
            row=0, column=2, sticky="w")

        self._bake_frame = tk.Frame(p, bg=self.pal["bg"])
        tk.Label(self._bake_frame, text="Target time (years)",
                 bg=self.pal["bg"], fg=self.pal["dim"],
                 font=FONT_LABEL).pack(fill="x", padx=14, pady=(4, 0))
        self.target_entry = tk.Entry(
            self._bake_frame, bg="#12141e", fg=self.pal["text"],
            insertbackground=self.pal["text"], relief="flat",
            highlightbackground=self.pal["border"],
            highlightthickness=1, font=FONT_LABEL)
        self.target_entry.insert(0, self.config.get("default_target", "10.0"))
        self.target_entry.pack(fill="x", padx=14, pady=(4, 4), ipady=4)
        self._on_rt_toggle()

    def _on_rt_toggle(self):
        if self.rt_var.get():
            self._bake_frame.pack_forget()
        else:
            self._bake_frame.pack(fill="x")

    def _build_panel_output(self, p: Panel):
        p.section("Output")
        p.label("File name  (no extension)")
        self.outfile_entry = p.entry(
            default=self.config.get("default_out_file", "frames"))
        p.separator()
        p.section("Result")
        self._result_box = tk.Text(
            p, bg="#0e1018", fg=self.pal["text"],
            relief="flat",
            highlightbackground=self.pal["border"],
            highlightthickness=1,
            font=("Courier New", 10),
            wrap="word", state="disabled",
        )
        self._result_box.pack(fill="both", expand=True, padx=14, pady=(4, 12))
        self._set_result("← Press Start to generate flags")

    def _set_result(self, text: str):
        self._result_box.configure(state="normal")
        self._result_box.delete("1.0", "end")
        self._result_box.insert("1.0", text)
        self._result_box.configure(state="disabled")

    # ── Нижняя полоса (прогресс + Start) ──────────────────
    def _build_bottom_bar(self):
        bar = tk.Frame(self, bg="#0e1018",
                       highlightbackground=self.pal["border"],
                       highlightthickness=1)
        bar.place(relx=0, rely=1.0, relwidth=1.0, height=54, anchor="sw")
        bar.grid_columnconfigure(0, weight=1)

        # статус
        self.status_lbl = tk.Label(bar, text="Ready to compute",
                                   bg="#0e1018", fg=self.pal["dim"],
                                   font=FONT_STATUS, anchor="w")
        self.status_lbl.grid(row=1, column=0, sticky="w", padx=16, pady=(0, 4))

        # прогресс-бар (Canvas)
        self._pb_canvas = tk.Canvas(bar, height=6, bg="#1e2030",
                                    highlightthickness=0)
        self._pb_canvas.grid(row=0, column=0, sticky="ew",
                             padx=16, pady=(8, 2))
        self._pb_fill = None
        self._pb_pct_lbl = tk.Label(bar, text="0%", bg="#0e1018",
                                    fg=self.pal["dim"],
                                    font=FONT_PROGRESS)
        self._pb_pct_lbl.grid(row=0, column=1, padx=(4, 8))

        # Start
        self._start_btn = RoundedButton(
            bar, text="Start", width=130, height=38,
            bg_color=self.pal["start"],
            bg_hover=self.pal["start_hover"],
            command=self._on_start,
        )
        self._start_btn.grid(row=0, column=2, rowspan=2,
                             padx=(0, 16), pady=6)

        self._pb_canvas.bind("<Configure>", self._redraw_pb)

    def _redraw_pb(self, _=None):
        self._pb_canvas.delete("all")
        w = self._pb_canvas.winfo_width()
        h = self._pb_canvas.winfo_height()
        r = h // 2
        # track
        self._pb_canvas.create_arc(0, 0, r*2, h, start=90, extent=180,
                                   fill="#2e3240", outline="")
        self._pb_canvas.create_arc(w - r*2, 0, w, h, start=270, extent=180,
                                   fill="#2e3240", outline="")
        self._pb_canvas.create_rectangle(r, 0, w - r, h,
                                         fill="#2e3240", outline="")
        # fill
        fw = int(w * self._progress)
        if fw > r * 2:
            self._pb_canvas.create_arc(0, 0, r*2, h, start=90, extent=180,
                                       fill=self.pal["start"], outline="")
            if fw < w - r:
                self._pb_canvas.create_rectangle(r, 0, fw, h,
                                                 fill=self.pal["start"], outline="")
            else:
                self._pb_canvas.create_rectangle(r, 0, w - r, h,
                                                 fill=self.pal["start"], outline="")
                self._pb_canvas.create_arc(w - r*2, 0, w, h, start=270,
                                           extent=180, fill=self.pal["start"],
                                           outline="")

    # ════════════════════════════════════════════════════════
    def _on_start(self):
        """Собрать флаги и показать в Output-панели."""
        src = self.source_var.get()
        src_idx = str(SOURCES.index(src))

        stdin: list[tuple[str, str]] = [
            (self.render_var.get(),  f"render_mode  →  {dict(RENDER_MODES).get(self.render_var.get(), '?')}... "
                                     f"(0=Cubes 1=Sphere 2=Points)"),
            (self.algo_var.get(),    f"algorithm    →  {dict(ALGORITHMS).get(self.algo_var.get(), '?')}"),
            (self.dt_entry.get(),    f"dt           →  {self.dt_entry.get()} years/step"),
            (src_idx,                f"source       →  {src}"),
        ]
        if src == SOURCES[0]:
            stdin.append((self.h5_var.get(), f"h5_file      →  {self.h5_var.get()}"))
        elif src == SOURCES[1]:
            stdin.append((self.txt_entry.get(), f"txt_path     →  {self.txt_entry.get()}"))
        else:
            stdin.append((self.rand_count.get(), f"body_count   →  {self.rand_count.get()}"))
            stdin.append((self.scenario_var.get(), f"scenario     →  {SCENARIOS[int(self.scenario_var.get())-1]}"))

        outfile = self.outfile_entry.get().strip() or "frames"
        is_rt   = self.rt_var.get()
        stdin += [
            (outfile,                f"output       →  data/{outfile}.h5"),
            (self.save_n_entry.get(), f"save_every   →  {self.save_n_entry.get()} steps"),
            ("1" if is_rt else "0",  f"realtime     →  {'yes' if is_rt else 'no (Bake)'}"),
        ]
        if not is_rt:
            stdin.append((self.target_entry.get(), f"target_time  →  {self.target_entry.get()} years"))

        lines = ["─── stdin sequence ───────────────────"]
        for val, desc in stdin:
            lines.append(f"  {val:<22}  # {desc}")
        lines += ["", "─── future CLI flags ─────────────────",
                  f"  ./Simulate --render {self.render_var.get()} \\",
                  f"    --algo {self.algo_var.get()} --dt {self.dt_entry.get()} \\",
                  f"    --source {src_idx} --output {outfile}"]
        self._set_result("\n".join(lines))
        self.status_lbl.configure(text="Flags generated  (launch not implemented yet)")
        self._progress = 0.0
        self._pb_pct_lbl.configure(text="0%")
        self._redraw_pb()

    # helpers
    def _list_h5(self) -> list[str]:
        files = []
        if os.path.isdir(self.data_dir):
            for f in sorted(os.listdir(self.data_dir)):
                if f.endswith(".h5"):
                    files.append(os.path.join("data", f))
        return files or ["(no .h5 files)"]

    def _browse_txt(self):
        path = filedialog.askopenfilename(
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
        if path:
            self.txt_entry.delete(0, "end")
            self.txt_entry.insert(0, path)
