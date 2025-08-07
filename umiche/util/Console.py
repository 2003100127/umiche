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

    def _stage(self, msg: str):
        """safe output (与 tqdm 兼容)"""
        from tqdm.auto import tqdm
        tqdm.write(msg)

    def _fmt_kv(self, **kw):
        return " | ".join(f"{k}={v}" for k, v in kw.items() if v is not None)

    # === 新增：对外暴露简洁版 k=v 拼接 ===
    def kv(self, **kw):
        return self._fmt_kv(**kw)

    def step(self, msg: str):
        if not self.verbose:
            return
        rid = Console._current_run_id.get()
        prefix = f"[{rid}] " if rid else ""
        # 用 _stage 保证与 tqdm 进度条兼容；若你更喜欢带时间戳，可改为 self.print
        self._stage(prefix + msg)

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