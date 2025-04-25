__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from functools import wraps
from umiche.plot.gadget.Draw import Draw as draw
from umiche.plot.gadget.Touch import Touch as touch


class Transmitter:

    def __init__(
            self,
            type,
            task=None,
    ):
        self.type = type
        self.task = task
        self.touch = touch()

    def __call__(self, deal):
        if self.type == "bar_individual":
            draw_canvas = self.bar_individual
        elif self.type == "hist_grouped":
            draw_canvas = self.hist_grouped
        elif self.type == "line":
            draw_canvas = self.line
        elif self.type == "line_scatter":
            draw_canvas = self.line_scatter
        if self.task == 'annotate_height':
            tinker = self.touch.bar_height_label
        else:
            tinker = None
        @wraps(deal)
        def canvas(self, *args, **kwargs):
            # print(kwargs)
            proc_f = draw_canvas(*args, **kwargs)
            # print(proc_f)
            if tinker:
                tinker(
                    canvas=proc_f['canvas'],
                    ax=proc_f['ax'],
                    format=proc_f['y_format'],
                    fontsize=10,
                )
            return deal(self, *args, **kwargs)
        return canvas

    @draw(type='bar_individual')
    def bar_individual(*args, **kwargs):
            return kwargs

    @draw(type='hist_grouped')
    def hist_grouped(*args, **kwargs):
        return kwargs

    @draw(type='line')
    def line(*args, **kwargs):
        return kwargs

    @draw(type='line_scatter')
    def line_scatter(*args, **kwargs):
        return kwargs