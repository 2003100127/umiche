__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

import networkx as nx
from umiche.network.CC import cc as gbfscc


class Cluster:


    def cc(
            self,
            graph_adj,
    ):
        """

        Parameters
        ----------
        graph_adj

        Returns
        -------

        """
        connected_components = list(gbfscc().deque(graph_adj))
        return {i: cc for i, cc in enumerate(connected_components)}

    def ccnx(
            self,
            edge_list,
    ):
        """

        Parameters
        ----------
        edge_list

        Returns
        -------

        """
        G = nx.Graph()
        for edge in edge_list:
            G.add_edge(edge[0], edge[1])
        return {i: G.subgraph(cc).nodes() for i, cc in enumerate(nx.connected_components(G))}


if __name__ == "__main__":
    p = Cluster()
    # print(p.cc({0: [1,], 1: [0, 2], 2: [1]}))
    # print(p.cc({0: []}))
    asds = p.cc(
        {0: [5, 9, 10, 25, 51, 56, 69, 79, 81, 91, 110],
         1: [4, 26, 98, 104, 122, 177, 184, 238, 251],
         2: [9, 23, 26, 48, 54, 56, 98, 223],
         3: [14, 18, 36, 42, 146, 165],
         4: [1, 17, 22, 26, 158, 237, 238],
         5: [0, 10, 36, 38, 51, 55, 59, 127, 146, 172],
         6: [34, 134, 149, 190, 222],
         7: [55, 65, 74, 82, 96, 165, 178],
         8: [10, 18, 20, 30, 36, 41, 81, 166, 182],
         9: [0, 2, 20, 33, 38, 47, 69, 81, 120, 159, 221, 223],
         10: [0, 5, 8, 20, 49, 51, 75, 115, 153, 166, 199],
         11: [13, 18, 21, 22, 37, 49, 96, 165, 166, 237],
         12: [185],
         13: [11, 16, 41, 75, 96, 122, 166, 184, 207],
         14: [3, 18, 21, 27, 31, 43, 50, 118, 119, 182],
         15: [70, 127, 205, 229, 240],
         16: [13, 20, 37, 41, 44, 47, 71, 75, 100, 140, 144, 241],
         17: [4, 25, 56, 111, 158, 177, 189, 229],
         18: [3, 8, 11, 14, 30, 37, 41, 49],
         19: [35, 137, 182, 196],
         20: [8, 9, 10, 16, 33, 34, 37, 38, 54, 136, 144, 166],
         21: [11, 14, 22, 31, 93, 119, 165, 178, 207],
         22: [4, 11, 21, 25, 65, 69, 122, 165, 216, 221],
         23: [2, 26, 56, 81, 97, 104, 105, 138, 158, 208, 209, 231],
         24: [121, 147, 150, 167, 199],
         25: [0, 17, 22, 49, 67, 79, 91, 119, 146, 189, 221],
         26: [1, 2, 4, 23, 56, 64, 69, 99, 102, 103, 126, 156, 204, 238, 244],
         27: [14, 28, 118, 158],
         28: [27, 31, 39, 171],
         29: [45, 173, 230],
         30: [8, 18, 41, 50, 80, 96, 115, 135, 144],
         31: [14, 21, 28, 33, 37, 44, 119, 129, 171, 212, 221],
         32: [48, 154, 179, 247],
         33: [9, 20, 31, 35, 38, 44, 51, 58, 93, 129, 182],
         34: [6, 20, 35, 54, 100, 134, 139, 160, 164, 223],
         35: [19, 33, 34, 139, 171, 196, 223],
         36: [3, 5, 8, 38, 42, 55, 81, 137, 182, 209],
         37: [11, 16, 18, 20, 31, 49, 144, 221],
         38: [5, 9, 20, 33, 36, 55, 71, 82, 139, 145, 206],
         39: [28, 44, 98, 133, 140, 155, 200],
         40: [117],
         41: [8, 13, 16, 18, 30, 42, 43, 68, 75, 90, 162, 220],
         42: [3, 36, 41, 43, 59, 71, 90, 132, 151, 232],
         43: [14, 41, 42, 44, 50, 53, 90, 182, 207],
         44: [16, 31, 33, 39, 43, 47, 53, 71, 129, 207],
         45: [29, 63, 64, 123, 175, 197, 204],
         46: [98, 107, 155, 174, 239, 247],
         47: [9, 16, 44, 71, 79, 90, 98, 113, 122, 159, 221, 248],
         48: [2, 32, 105, 149, 154, 170, 204, 222, 247],
         49: [10, 11, 18, 25, 37, 75, 87, 115, 119, 146],
         50: [14, 30, 43, 80, 129, 178, 182],
         51: [0, 5, 10, 33, 53, 60, 93, 94, 119, 182, 196],
         52: [83, 110, 152, 170],
         53: [43, 44, 51, 59, 75, 76, 79, 119, 133, 207],
         54: [2, 20, 34, 140, 149, 199, 244],
         55: [5, 7, 36, 38, 63, 64, 69, 93, 95, 165, 166, 227],
         56: [0, 2, 17, 23, 26, 73, 94, 127, 147, 170, 177, 199, 229],
         57: [215, 217, 243],
         58: [33, 145],
         59: [5, 42, 53, 71, 75, 79, 146, 205],
         60: [51, 121, 224],
         61: [63, 66, 69, 83, 109, 110, 194, 204, 216],
         62: [],
         63: [45, 55, 61, 66, 145, 173],
         64: [26, 45, 55, 95, 102, 127, 209, 227, 244],
         65: [7, 22, 69, 80, 91, 96, 122, 159, 178, 238],
         66: [61, 63, 123, 166, 225],
         67: [25, 108, 110, 111, 216],
         68: [41, 100, 143, 162, 220],
         69: [0, 9, 22, 26, 55, 61, 65, 81, 93, 103, 122, 126, 166, 194],
         70: [15, 82, 179, 200, 240],
         71: [16, 38, 42, 44, 47, 59, 82, 114, 157, 176, 200],
         72: [100, 134, 193, 241],
         73: [56, 97, 147, 156, 170, 192, 215],
         74: [7, 195, 227],
         75: [10, 13, 16, 41, 49, 53, 59, 79, 115, 191],
         76: [53, 133, 161, 196, 233],
         77: [88, 106],
         78: [],
         79: [0, 25, 47, 53, 59, 75, 90, 91, 122, 161, 177],
         80: [30, 50, 65, 81, 90, 91, 101, 138, 159, 187],
         81: [0, 8, 9, 23, 36, 69, 80, 90, 182, 231],
         82: [7, 38, 70, 71, 128, 129, 144, 159],
         83: [52, 61, 92, 103, 204],
         84: [232],
         85: [89, 198, 236, 251],
         86: [],
         87: [49],
         88: [77, 105, 108, 111, 154, 158],
         89: [85, 143, 174, 198, 246],
         90: [41, 42, 43, 47, 79, 80, 81, 104, 122, 143],
         91: [0, 25, 65, 79, 80, 115, 159, 229],
         92: [83, 99, 103, 194],
         93: [21, 33, 51, 55, 69, 102, 166, 178, 182, 207],
         94: [51, 56, 102, 121, 127, 133, 192, 196, 199],
         95: [55, 64, 126, 137, 139, 164, 168, 214, 227],
         96: [7, 11, 13, 30, 65, 115, 144, 166, 178, 195, 226, 245],
         97: [23, 73, 105, 156, 163, 208],
         98: [1, 2, 39, 46, 47, 104, 113, 140, 177, 200, 247, 248],
         99: [26, 92, 147, 156, 180, 194, 204, 208],
         100: [16, 34, 68, 72, 114, 140, 160, 248],
         101: [80, 187],
         102: [26, 64, 93, 94, 175, 180, 228, 244],
         103: [26, 69, 83, 92, 126, 227, 231],
         104: [1, 23, 90, 98, 132, 138, 143, 158, 162, 174, 177],
         105: [23, 48, 88, 97, 170, 204, 208],
         106: [77, 162, 181],
         107: [46, 140, 155, 239],
         108: [67, 88, 187, 216],
         109: [61, 122, 173, 216, 236],
         110: [0, 52, 61, 67, 152, 170],
         111: [17, 67, 88, 117, 141, 154, 170],
         112: [210],
         113: [47, 98, 202, 248],
         114: [71, 100, 139, 200, 219, 248],
         115: [10, 30, 49, 75, 91, 96, 144],
         116: [145, 183, 186, 190, 206],
         117: [40, 111, 131],
         118: [14, 27],
         119: [14, 21, 25, 31, 49, 51, 53, 146],
         120: [9, 194],
         121: [24, 60, 94, 147, 167, 180, 192, 224],
         122: [1, 13, 22, 47, 65, 69, 79, 90, 109, 207, 236, 251],
         123: [45, 66, 149, 150, 175, 204, 244],
         124: [],
         125: [154, 201],
         126: [26, 69, 95, 103, 164, 203, 223, 251],
         127: [5, 15, 56, 64, 94, 167, 172, 197, 199, 205, 209],
         128: [82, 145, 176, 179],
         129: [31, 33, 44, 50, 82, 144, 159, 178],
         130: [213],
         131: [117],
         132: [42, 104, 151, 162, 200, 205, 209, 240],
         133: [39, 53, 76, 94, 177, 191, 205],
         134: [6, 34, 72, 136],
         135: [30, 138, 162, 240],
         136: [20, 134, 153, 241],
         137: [19, 36, 95, 139, 209, 242],
         138: [23, 80, 104, 135, 158, 229, 238, 240],
         139: [34, 35, 38, 95, 114, 137, 206, 223],
         140: [16, 39, 54, 98, 100, 107, 162, 184, 191, 200],
         141: [111, 170, 177, 247],
         142: [183, 186, 210],
         143: [68, 89, 90, 104, 161, 246, 248, 251],
         144: [16, 20, 30, 37, 82, 96, 115, 129, 159, 160],
         145: [38, 58, 63, 116, 128, 176],
         146: [3, 5, 25, 49, 59, 119, 165, 211],
         147: [24, 56, 73, 99, 121, 167, 170, 208],
         148: [211],
         149: [6, 48, 54, 123, 150, 190],
         150: [24, 123, 149, 170, 197, 199],
         151: [42, 132, 220],
         152: [52, 110, 170, 215, 222, 243],
         153: [10, 136],
         154: [32, 48, 88, 111, 125, 247],
         155: [39, 46, 107, 213, 239, 250],
         156: [26, 73, 97, 99, 203, 204],
         157: [71, 176, 232, 241],
         158: [4, 17, 23, 27, 88, 104, 138],
         159: [9, 47, 65, 80, 82, 91, 129, 144, 221],
         160: [34, 100, 144],
         161: [76, 79, 143, 177, 198, 248, 251],
         162: [41, 68, 104, 106, 132, 135, 140, 184, 191, 220],
         163: [97],
         164: [34, 95, 126, 166, 244],
         165: [3, 7, 11, 21, 22, 55, 146, 168],
         166: [8, 10, 11, 13, 20, 55, 66, 69, 93, 96, 164, 225, 244],
         167: [24, 121, 127, 147, 197, 217],
         168: [95, 165],
         169: [175],
         170: [48, 52, 56, 73, 105, 110, 111, 141, 147, 150, 152, 197, 204],
         171: [28, 31, 35, 188],
         172: [5, 127, 206, 227],
         173: [29, 63, 109, 176],
         174: [46, 89, 104, 208],
         175: [45, 102, 123, 169, 180, 204],
         176: [71, 128, 145, 157, 173, 219, 230],
         177: [1, 17, 56, 79, 98, 104, 133, 141, 161, 191, 205, 229],
         178: [7, 21, 50, 65, 93, 96, 129, 207, 228],
         179: [32, 70, 128, 230],
         180: [99, 102, 121, 175],
         181: [106],
         182: [8, 14, 19, 33, 36, 43, 50, 51, 81, 93],
         183: [116, 142, 210],
         184: [1, 13, 140, 162, 191, 237, 244],
         185: [12],
         186: [116, 142, 206],
         187: [80, 101, 108],
         188: [171, 202, 221, 223, 248],
         189: [17, 25, 202],
         190: [6, 116, 149],
         191: [75, 133, 140, 162, 177, 184, 199, 205],
         192: [73, 94, 121],
         193: [72],
         194: [61, 69, 92, 99, 120, 225, 236],
         195: [74, 96],
         196: [19, 35, 51, 76, 94, 224, 233],
         197: [45, 127, 150, 167, 170, 243],
         198: [85, 89, 161],
         199: [10, 24, 54, 56, 94, 127, 150, 191, 244],
         200: [39, 70, 71, 98, 114, 132, 140, 205, 230, 239],
         201: [125],
         202: [113, 188, 189, 221],
         203: [126, 156, 215, 218],
         204: [26, 45, 48, 61, 83, 99, 105, 123, 156, 170, 175],
         205: [15, 59, 127, 132, 133, 177, 191, 200],
         206: [38, 116, 139, 172, 186, 227],
         207: [13, 21, 43, 44, 53, 93, 122, 178],
         208: [23, 97, 99, 105, 147, 174],
         209: [23, 36, 64, 127, 132, 137, 240],
         210: [112, 142, 183, 242],
         211: [146, 148],
         212: [31],
         213: [130, 155, 250],
         214: [95, 217, 242],
         215: [57, 73, 152, 203],
         216: [22, 61, 67, 108, 109],
         217: [57, 167, 214, 224, 242, 243],
         218: [203],
         219: [114, 176, 230],
         220: [41, 68, 151, 162],
         221: [9, 22, 25, 31, 37, 47, 159, 188, 202],
         222: [6, 48, 152, 223],
         223: [2, 9, 34, 35, 126, 139, 188, 222, 248],
         224: [60, 121, 196, 217],
         225: [66, 166, 194, 245],
         226: [96, 245],
         227: [55, 64, 74, 95, 103, 172, 206],
         228: [102, 178, 238],
         229: [15, 17, 56, 91, 138, 177, 238],
         230: [29, 176, 179, 200, 219, 239, 247],
         231: [23, 81, 103],
         232: [42, 84, 157],
         233: [76, 196],
         234: [],
         235: [],
         236: [85, 109, 122, 194],
         237: [4, 11, 184, 244],
         238: [1, 4, 26, 65, 138, 228, 229],
         239: [46, 107, 155, 200, 230],
         240: [15, 70, 132, 135, 138, 209],
         241: [16, 72, 136, 157],
         242: [137, 210, 214, 217],
         243: [57, 152, 197, 217],
         244: [26, 54, 64, 102, 123, 164, 166, 184, 199, 237],
         245: [96, 225, 226],
         246: [89, 143],
         247: [32, 46, 48, 98, 141, 154, 230],
         248: [47, 98, 100, 113, 114, 143, 161, 188, 223, 251],
         249: [],
         250: [155, 213],
         251: [1, 85, 122, 126, 143, 161, 248],
         })
    print(asds)
    print(p.ccnx([[0, 5], [0, 9], [0, 10], [0, 25], [0, 51], [0, 56], [0, 69], [0, 79], [0, 81], [0, 91], [0, 110], [1, 4], [1, 26], [1, 98], [1, 104], [1, 122], [1, 177], [1, 184], [1, 238], [1, 251], [2, 9], [2, 23], [2, 26], [2, 48], [2, 54], [2, 56], [2, 98], [2, 223], [3, 14], [3, 18], [3, 36], [3, 42], [3, 146], [3, 165], [4, 17], [4, 22], [4, 26], [4, 158], [4, 237], [4, 238], [5, 10], [5, 36], [5, 38], [5, 51], [5, 55], [5, 59], [5, 127], [5, 146], [5, 172], [6, 34], [6, 134], [6, 149], [6, 190], [6, 222], [7, 55], [7, 65], [7, 74], [7, 82], [7, 96], [7, 165], [7, 178], [8, 10], [8, 18], [8, 20], [8, 30], [8, 36], [8, 41], [8, 81], [8, 166], [8, 182], [9, 20], [9, 33], [9, 38], [9, 47], [9, 69], [9, 81], [9, 120], [9, 159], [9, 221], [9, 223], [10, 20], [10, 49], [10, 51], [10, 75], [10, 115], [10, 153], [10, 166], [10, 199], [11, 13], [11, 18], [11, 21], [11, 22], [11, 37], [11, 49], [11, 96], [11, 165], [11, 166], [11, 237], [12, 185], [13, 16], [13, 41], [13, 75], [13, 96], [13, 122], [13, 166], [13, 184], [13, 207], [14, 18], [14, 21], [14, 27], [14, 31], [14, 43], [14, 50], [14, 118], [14, 119], [14, 182], [15, 70], [15, 127], [15, 205], [15, 229], [15, 240], [16, 20], [16, 37], [16, 41], [16, 44], [16, 47], [16, 71], [16, 75], [16, 100], [16, 140], [16, 144], [16, 241], [17, 25], [17, 56], [17, 111], [17, 158], [17, 177], [17, 189], [17, 229], [18, 30], [18, 37], [18, 41], [18, 49], [19, 35], [19, 137], [19, 182], [19, 196], [20, 33], [20, 34], [20, 37], [20, 38], [20, 54], [20, 136], [20, 144], [20, 166], [21, 22], [21, 31], [21, 93], [21, 119], [21, 165], [21, 178], [21, 207], [22, 25], [22, 65], [22, 69], [22, 122], [22, 165], [22, 216], [22, 221], [23, 26], [23, 56], [23, 81], [23, 97], [23, 104], [23, 105], [23, 138], [23, 158], [23, 208], [23, 209], [23, 231], [24, 121], [24, 147], [24, 150], [24, 167], [24, 199], [25, 49], [25, 67], [25, 79], [25, 91], [25, 119], [25, 146], [25, 189], [25, 221], [26, 56], [26, 64], [26, 69], [26, 99], [26, 102], [26, 103], [26, 126], [26, 156], [26, 204], [26, 238], [26, 244], [27, 28], [27, 118], [27, 158], [28, 31], [28, 39], [28, 171], [29, 45], [29, 173], [29, 230], [30, 41], [30, 50], [30, 80], [30, 96], [30, 115], [30, 135], [30, 144], [31, 33], [31, 37], [31, 44], [31, 119], [31, 129], [31, 171], [31, 212], [31, 221], [32, 48], [32, 154], [32, 179], [32, 247], [33, 35], [33, 38], [33, 44], [33, 51], [33, 58], [33, 93], [33, 129], [33, 182], [34, 35], [34, 54], [34, 100], [34, 134], [34, 139], [34, 160], [34, 164], [34, 223], [35, 139], [35, 171], [35, 196], [35, 223], [36, 38], [36, 42], [36, 55], [36, 81], [36, 137], [36, 182], [36, 209], [37, 49], [37, 144], [37, 221], [38, 55], [38, 71], [38, 82], [38, 139], [38, 145], [38, 206], [39, 44], [39, 98], [39, 133], [39, 140], [39, 155], [39, 200], [40, 117], [41, 42], [41, 43], [41, 68], [41, 75], [41, 90], [41, 162], [41, 220], [42, 43], [42, 59], [42, 71], [42, 90], [42, 132], [42, 151], [42, 232], [43, 44], [43, 50], [43, 53], [43, 90], [43, 182], [43, 207], [44, 47], [44, 53], [44, 71], [44, 129], [44, 207], [45, 63], [45, 64], [45, 123], [45, 175], [45, 197], [45, 204], [46, 98], [46, 107], [46, 155], [46, 174], [46, 239], [46, 247], [47, 71], [47, 79], [47, 90], [47, 98], [47, 113], [47, 122], [47, 159], [47, 221], [47, 248], [48, 105], [48, 149], [48, 154], [48, 170], [48, 204], [48, 222], [48, 247], [49, 75], [49, 87], [49, 115], [49, 119], [49, 146], [50, 80], [50, 129], [50, 178], [50, 182], [51, 53], [51, 60], [51, 93], [51, 94], [51, 119], [51, 182], [51, 196], [52, 83], [52, 110], [52, 152], [52, 170], [53, 59], [53, 75], [53, 76], [53, 79], [53, 119], [53, 133], [53, 207], [54, 140], [54, 149], [54, 199], [54, 244], [55, 63], [55, 64], [55, 69], [55, 93], [55, 95], [55, 165], [55, 166], [55, 227], [56, 73], [56, 94], [56, 127], [56, 147], [56, 170], [56, 177], [56, 199], [56, 229], [57, 215], [57, 217], [57, 243], [58, 145], [59, 71], [59, 75], [59, 79], [59, 146], [59, 205], [60, 121], [60, 224], [61, 63], [61, 66], [61, 69], [61, 83], [61, 109], [61, 110], [61, 194], [61, 204], [61, 216], [63, 66], [63, 145], [63, 173], [64, 95], [64, 102], [64, 127], [64, 209], [64, 227], [64, 244], [65, 69], [65, 80], [65, 91], [65, 96], [65, 122], [65, 159], [65, 178], [65, 238], [66, 123], [66, 166], [66, 225], [67, 108], [67, 110], [67, 111], [67, 216], [68, 100], [68, 143], [68, 162], [68, 220], [69, 81], [69, 93], [69, 103], [69, 122], [69, 126], [69, 166], [69, 194], [70, 82], [70, 179], [70, 200], [70, 240], [71, 82], [71, 114], [71, 157], [71, 176], [71, 200], [72, 100], [72, 134], [72, 193], [72, 241], [73, 97], [73, 147], [73, 156], [73, 170], [73, 192], [73, 215], [74, 195], [74, 227], [75, 79], [75, 115], [75, 191], [76, 133], [76, 161], [76, 196], [76, 233], [77, 88], [77, 106], [79, 90], [79, 91], [79, 122], [79, 161], [79, 177], [80, 81], [80, 90], [80, 91], [80, 101], [80, 138], [80, 159], [80, 187], [81, 90], [81, 182], [81, 231], [82, 128], [82, 129], [82, 144], [82, 159], [83, 92], [83, 103], [83, 204], [84, 232], [85, 89], [85, 198], [85, 236], [85, 251], [88, 105], [88, 108], [88, 111], [88, 154], [88, 158], [89, 143], [89, 174], [89, 198], [89, 246], [90, 104], [90, 122], [90, 143], [91, 115], [91, 159], [91, 229], [92, 99], [92, 103], [92, 194], [93, 102], [93, 166], [93, 178], [93, 182], [93, 207], [94, 102], [94, 121], [94, 127], [94, 133], [94, 192], [94, 196], [94, 199], [95, 126], [95, 137], [95, 139], [95, 164], [95, 168], [95, 214], [95, 227], [96, 115], [96, 144], [96, 166], [96, 178], [96, 195], [96, 226], [96, 245], [97, 105], [97, 156], [97, 163], [97, 208], [98, 104], [98, 113], [98, 140], [98, 177], [98, 200], [98, 247], [98, 248], [99, 147], [99, 156], [99, 180], [99, 194], [99, 204], [99, 208], [100, 114], [100, 140], [100, 160], [100, 248], [101, 187], [102, 175], [102, 180], [102, 228], [102, 244], [103, 126], [103, 227], [103, 231], [104, 132], [104, 138], [104, 143], [104, 158], [104, 162], [104, 174], [104, 177], [105, 170], [105, 204], [105, 208], [106, 162], [106, 181], [107, 140], [107, 155], [107, 239], [108, 187], [108, 216], [109, 122], [109, 173], [109, 216], [109, 236], [110, 152], [110, 170], [111, 117], [111, 141], [111, 154], [111, 170], [112, 210], [113, 202], [113, 248], [114, 139], [114, 200], [114, 219], [114, 248], [115, 144], [116, 145], [116, 183], [116, 186], [116, 190], [116, 206], [117, 131], [119, 146], [120, 194], [121, 147], [121, 167], [121, 180], [121, 192], [121, 224], [122, 207], [122, 236], [122, 251], [123, 149], [123, 150], [123, 175], [123, 204], [123, 244], [125, 154], [125, 201], [126, 164], [126, 203], [126, 223], [126, 251], [127, 167], [127, 172], [127, 197], [127, 199], [127, 205], [127, 209], [128, 145], [128, 176], [128, 179], [129, 144], [129, 159], [129, 178], [130, 213], [132, 151], [132, 162], [132, 200], [132, 205], [132, 209], [132, 240], [133, 177], [133, 191], [133, 205], [134, 136], [135, 138], [135, 162], [135, 240], [136, 153], [136, 241], [137, 139], [137, 209], [137, 242], [138, 158], [138, 229], [138, 238], [138, 240], [139, 206], [139, 223], [140, 162], [140, 184], [140, 191], [140, 200], [141, 170], [141, 177], [141, 247], [142, 183], [142, 186], [142, 210], [143, 161], [143, 246], [143, 248], [143, 251], [144, 159], [144, 160], [145, 176], [146, 165], [146, 211], [147, 167], [147, 170], [147, 208], [148, 211], [149, 150], [149, 190], [150, 170], [150, 197], [150, 199], [151, 220], [152, 170], [152, 215], [152, 222], [152, 243], [154, 247], [155, 213], [155, 239], [155, 250], [156, 203], [156, 204], [157, 176], [157, 232], [157, 241], [159, 221], [161, 177], [161, 198], [161, 248], [161, 251], [162, 184], [162, 191], [162, 220], [164, 166], [164, 244], [165, 168], [166, 225], [166, 244], [167, 197], [167, 217], [169, 175], [170, 197], [170, 204], [171, 188], [172, 206], [172, 227], [173, 176], [174, 208], [175, 180], [175, 204], [176, 219], [176, 230], [177, 191], [177, 205], [177, 229], [178, 207], [178, 228], [179, 230], [183, 210], [184, 191], [184, 237], [184, 244], [186, 206], [188, 202], [188, 221], [188, 223], [188, 248], [189, 202], [191, 199], [191, 205], [194, 225], [194, 236], [196, 224], [196, 233], [197, 243], [199, 244], [200, 205], [200, 230], [200, 239], [202, 221], [203, 215], [203, 218], [206, 227], [209, 240], [210, 242], [213, 250], [214, 217], [214, 242], [217, 224], [217, 242], [217, 243], [219, 230], [222, 223], [223, 248], [225, 245], [226, 245], [228, 238], [229, 238], [230, 239], [230, 247], [237, 244], [248, 251]]))

