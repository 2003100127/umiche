__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from datetime import datetime

import time, uuid, inspect, functools, contextvars


class Console:

    # run-id for step as prefix
    _current_run_id = contextvars.ContextVar("vignette_run_id", default=None)
    _current_func   = contextvars.ContextVar("vignette_func",    default=None)

    def __init__(
            self,
            placeholder='logger: ',
            verbose=False,
    ):
        self._verbose = verbose
        self.placeholder = placeholder

    @property
    def verbose(self, ):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = value

    def print(self, content):
        if self._verbose:
            now = datetime.now()
            dt_format = now.strftime("%d/%m/%Y %H:%M:%S ")
            print(dt_format + self.placeholder + str(content))

    def check(self, content):
        if self._verbose:
            print(content)

    def df_column_summary(
            self,
            df,
            title: str = "DataFrame Columns",
    ):
        """
        Display DataFrame column names in a colorful table with index numbers.

        Parameters
        ----------
        df
            (pd.DataFrame): The DataFrame whose columns to display.
        title
            (str): Title of the table.

        Returns
        -------

        """
        from rich.console import Console as RichConsole
        from rich.table import Table
        import pandas as pd
        import numpy as np
        import json

        # ---- helpers (do not change the external naming, only add robust tools within the function) ----
        def _preview_value(x, maxlen=30):
            """Provide refined previews for complex objects to avoid flooding the console with lengthy content."""
            if x is None or (isinstance(x, float) and pd.isna(x)):
                s = "<NaN>"
            elif isinstance(x, dict):
                try:
                    keys = list(x.keys())
                except Exception:
                    keys = []
                s = f"<dict len={len(x)} keys={keys[:3]}>"
            elif isinstance(x, (list, tuple, set)):
                try:
                    sample = list(x)[:3]
                    tname = type(x).__name__
                    s = f"<{tname} len={len(x)} sample={sample}>"
                except Exception:
                    s = f"<{type(x).__name__}>"
            elif isinstance(x, np.ndarray):
                s = f"<ndarray shape={x.shape} dtype={x.dtype}>"
            else:
                s = str(x)
            return s if len(s) <= maxlen else s[:maxlen - 3] + "..."

        def _to_hashable_for_unique(x):
            """Convert unhashable objects to stable hashable representations for safe nunique computation."""
            if x is None or (isinstance(x, float) and pd.isna(x)):
                return x
            if isinstance(x, np.ndarray):
                return tuple(x.tolist())
            if isinstance(x, (list, tuple)):
                return tuple(_to_hashable_for_unique(i) for i in x)
            if isinstance(x, set):
                return tuple(sorted(_to_hashable_for_unique(i) for i in x))
            if isinstance(x, dict):
                try:
                    return json.dumps(x, sort_keys=True, ensure_ascii=False, default=str)
                except Exception:
                    x2 = {str(_to_hashable_for_unique(k)): _to_hashable_for_unique(v) for k, v in x.items()}
                    return json.dumps(x2, sort_keys=True, ensure_ascii=False, default=str)
            try:
                hash(x)
                return x
            except Exception:
                return str(x)

        def _safe_nunique(series: pd.Series) -> int:
            """Safely compute the number of unique values for columns containing dict/list/ndarray etc."""
            try:
                return int(series.nunique(dropna=True))
            except TypeError:
                return int(series.map(_to_hashable_for_unique).nunique(dropna=True))

        rich_cons = RichConsole()
        table = Table(title=title + ': ' + str(df.shape[0]) + ' rows, ' + str(df.shape[1]) + ' columns', show_header=True, header_style="bold magenta")
        table.add_column("Index", justify="right", style="cyan", no_wrap=True)
        table.add_column("Column Name", style="green")
        table.add_column("First Value", style="yellow")
        table.add_column("Dtype", style="blue")
        table.add_column("Has NaN", style="red")
        table.add_column("Unique Count", style="magenta", justify="right", no_wrap=True)
        table.add_column("% NaN", style="red", justify="right")
        # table.add_column("Top3 Uniq Values", style="white")

        for i, col in enumerate(df.columns, start=1):
            # First value
            if not df.empty:
                first_val_raw = df[col].iloc[0]
                # Preview complex objects; preserve the original semantics of "display <NaN> if the value is NaN"
                first_val = _preview_value(first_val_raw, maxlen=30)
            else:
                first_val = "<empty>"

            # Additional truncation (consistent with the original logic to ensure it doesn't end up too long)
            if len(first_val) > 30:
                first_val = first_val[:27] + "..."

            # Data type
            dtype = str(df[col].dtype)

            # NaN check
            nan_count = int(df[col].isna().sum())
            has_nan = "Yes" if nan_count > 0 else "No"
            nan_pct = f"{(nan_count / len(df) * 100):.1f}%" if len(df) > 0 else "0.0%"

            # Unique value count
            unique_count = str(_safe_nunique(df[col]))

            # Top3 Uniq Values
            # examples = df[col].dropna().unique()[:3]
            # examples_str = ", ".join(map(str, examples)) if len(examples) > 0 else "<No Data>"
            # if len(examples_str) > 40:
            #     examples_str = examples_str[:37] + "..."

            table.add_row(
                str(i),
                col,
                first_val,
                dtype,
                has_nan,
                unique_count,
                nan_pct,
                # examples_str,
            )

        rich_cons.print(table)

    def _tqdm(
            self,
            iterable,
            desc: str = "",
            total: int | None = None,
            position: int = 0,
            leave: bool = False,
            unit: str = "unit",
            dynamic_ncols: bool = True,
    ):
        from tqdm.auto import tqdm
        disable = not bool(getattr(self, "verbose", True))

        return tqdm(
            iterable,
            total=total,
            desc=desc,
            position=position,
            leave=leave,
            dynamic_ncols=dynamic_ncols,
            unit=unit,
            disable=disable,
        )

    @staticmethod
    def vignette(label: str | None = None):
        """
        装饰“实例方法”（第一个参数为 self）。
        会优先使用被装饰实例上的 self.console；若不存在则临时创建一个 Console(verbose=True)。
        """
        def deco(func):
            @functools.wraps(func)
            def wrapper(self, *args, **kwargs):
                console = getattr(self, "console", None)
                if console is None:
                    console = Console(verbose=True)  # fallback：确保可见

                rid = uuid.uuid4().hex[:6].upper()
                Console._current_run_id.set(rid)
                Console._current_func.set(f"{self.__class__.__name__}.{func.__name__}")

                if console.verbose:
                    name = label or func.__name__
                    cls  = self.__class__.__name__
                    file = inspect.getsourcefile(func) or "<?>"
                    try:
                        line = inspect.getsourcelines(func)[1]
                    except Exception:
                        line = -1
                    hdr = f"▶ {cls}.{name} | id={rid} | method={getattr(self, 'clustering_method', None)}"
                    bar = "═" * max(10, min(100, len(hdr) + 8))
                    console.print(f"\n{bar}\n{hdr}\n↳ @ {file}:{line}\n{bar}")

                t0 = time.perf_counter()
                try:
                    return func(self, *args, **kwargs)
                finally:
                    if console.verbose:
                        dt = (time.perf_counter() - t0) * 1000
                        console.print(f"✔ DONE {self.__class__.__name__}.{func.__name__} in {dt:.1f} ms")
                    # clean context
                    Console._current_run_id.set(None)
                    Console._current_func.set(None)
            return wrapper
        return deco

    # === global staticmethod  ===
    @staticmethod
    def vignette_global(console=None, label: str | None = None):
        """
        used when class instance or self is not used, taking "console" as input whenever needed
        """
        def deco(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                cons = console or Console(verbose=True)
                rid  = uuid.uuid4().hex[:6].upper()
                name = label or func.__name__
                file = inspect.getsourcefile(func) or "<?>"
                try:
                    line = inspect.getsourcelines(func)[1]
                except Exception:
                    line = -1
                if cons.verbose:
                    bar = "═" * max(10, min(100, len(name) + 20))
                    cons.print(f"\n{bar}\n▶ {name} | id={rid}\n↳ @ {file}:{line}\n{bar}")
                t0 = time.perf_counter()
                try:
                    return func(*args, **kwargs)
                finally:
                    if cons.verbose:
                        dt = (time.perf_counter() - t0) * 1000
                        cons.print(f"✔ DONE {name} in {dt:.1f} ms")
            return wrapper
        return deco