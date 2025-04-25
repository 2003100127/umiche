__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


import numpy as np
import matplotlib.colors as mcolors


class Element:

    def __init__(self, ):
        self.color_list = self.color()

    def color(
            self,
            which='ccs4',
            is_random=False,
    ):
        """
        References
        ----------
            https://matplotlib.org/stable/gallery/color/named_colors.html

        It includes the color of the element
            'firebrick',
            'royalblue',
            'springgreen',
            'burlywood',
            'oldlace',
            'lightgoldenrodyellow',
            'grey',
            'cyan',
            'crimson',
            'mediumvioletred',
            'maroon',
            'mediumturquoise',
            'teal',
            'azure',
            'palevioletred',
            'mediumslateblue',
            'olivedrab',
            'darkred',
            'dimgrey',
            'lightgreen',
            'blueviolet',
            'chartreuse',
            'whitesmoke',
            'mediumspringgreen',
            'yellowgreen',
            'green',
            'darkorchid',
            'lawngreen',
            'black',
            'darkgoldenrod',
            'moccasin',
            'seashell',
            'orange',
            'slateblue',
            'hotpink',
            'tan',
            'darkolivegreen',
            'mediumorchid',
            'darkgreen',
            'navajowhite',
            'khaki',
            'paleturquoise',
            'darkgray',
            'darkorange',
            'salmon',
            'floralwhite',
            'lightgray',
            'darkblue',
            'lightslategray',
            'darkslateblue',
            'goldenrod',
            'lightcyan',
            'papayawhip',
            'wheat',
            'lightslategrey',
            'darksalmon',
            'skyblue',
            'violet',
            'coral',
            'mediumseagreen',
            'darkmagenta',
            'palegreen',
            'magenta',
            'lavenderblush',
            'darkkhaki',
            'darkturquoise',
            'lightsalmon',
            'blanchedalmond',
            'deepskyblue',
            'mediumaquamarine',
            'mintcream',
            'lightyellow',
            'linen',
            'peachpuff',
            'aliceblue',
            'gold',
            'ghostwhite',
            'silver',
            'bisque',
            'aquamarine',
            'lightblue',
            'lavender',
            'antiquewhite',
            'dimgray',
            'lightskyblue',
            'pink',
            'beige',
            'indianred',
            'rosybrown',
            'rebeccapurple',
            'lightgrey',
            'powderblue',
            'chocolate',
            'darkgrey',
            'purple',
            'darkslategrey',
            'cornsilk',
            'turquoise',
            'blue',
            'greenyellow',
            'red',
            'lime',
            'ivory',
            'royalblue',
            'palegoldenrod',
            'lightsteelblue',
            'slategray',
            'gray',
            'brown',
            'mistyrose',
            'cornflowerblue',
            'tomato',
            'indigo',
            'snow',
            'darkslategray',
            'lightpink',
            'darkviolet',
            'limegreen',
            'honeydew',
            'yellow',
            'sienna',
            'mediumblue',
            'fuchsia',
            'forestgreen',
            'orchid',
            'plum',
            'lemonchiffon',
            'gainsboro',
            'sandybrown',
            'slategrey',
            'mediumpurple',
            'saddlebrown',
            'lightcoral',
            'midnightblue',
            'navy',
            'thistle',
            'dodgerblue',
            'lightseagreen',
            'darkseagreen',
            'cadetblue',
            'seagreen',
            'darkcyan',
            'steelblue',
            'white',
            'peru',
            'deeppink',
            'olive',
            'firebrick',
            'aqua',
            'orangered',
        Returns
        -------

        """
        if which == 'base':
            colors = mcolors.BASE_COLORS
        elif which == 'ccs4':
            colors = mcolors.CSS4_COLORS
        elif which == 'tableau':
            colors = mcolors.TABLEAU_COLORS
        elif which == 'xkcd':
            colors = mcolors.XKCD_COLORS
        else:
            colors = mcolors.CSS4_COLORS
        # return np.random.permutation([*colors.keys()]).tolist()
        return [*colors.keys()]