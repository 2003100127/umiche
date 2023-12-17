__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "MIT"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"
__lab__ = "Cribbslab"

from collections import deque


class ConnectedComponent:

    def deque(
            self,
            graph,
    ):
        """

        Parameters
        ----------
        graph

        Returns
        -------

        """
        visited = set()
        for root, nbrs in graph.items():
            if root not in visited:
                visited.add(root)
                component = []
                queue = deque([root])
                # print('-> root {} has not been visited'.format(root))
                # print('===> a queue built by root {} is {}'.format(root, queue))
                while queue:
                    # print('======> a queue built by each root node {}'.format(queue))
                    node = queue.popleft()
                    # print('======> node: {}'.format(node))
                    component.append(node)
                    for nbr in graph[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
                # print('======> visited nodes {}'.format(visited))
                yield component
            else:
                continue
                # print('-> root {} has been visited'.format(root))

    def set(
            self,
            graph,
    ):
        """
        Examples
        --------
        graph_adj = {
            0: [4],
            1: [2, 5],
            2: [1, 3, 8],
            3: [2, 14],
            4: [0, 9],
            5: [1],
            6: [10, 12],
            7: [13],
            8: [2, 14],
            9: [4, 10, 15],
            10: [6, 9],
            11: [17],
            12: [6],
            13: [7, 19, 20],
            14: [3, 8],
            15: [9, 21],
            16: [22],
            17: [11, 18],
            18: [17, 19],
            19: [13, 18, 26],
            20: [13, 26],
            21: [15, 27],
            22: [16, 23],
            23: [22, 24, 28],
            24: [23, 25, 29],
            25: [24],
            26: [19, 20, 30, 31],
            27: [21],
            28: [23, 29],
            29: [24, 28],
            30: [26, 31],
            31: [26, 30],
        }

        ..  @ex:
            seen = set()
            def component(node):
                nodes = set([node])
                while nodes:
                    node = nodes.pop()
                    seen.add(node)
                    nodes |= neighbors[node] - seen
                    yield node
            for node in neighbors:
                if node not in seen:
                    yield component(node)

        Parameters
        ----------
        graph

        Returns
        -------

        """
        visited = set()
        components = []
        for root, nbrs in graph.items():
            if root not in visited:
                visited.add(root)
                component = []
                queue = [root]
                print('-> root {} has not been visited'.format(root))
                print('===> a queue built by root {} is {}'.format(root, queue))
                while queue:
                    print('======> a queue built by each root node {}'.format(queue))
                    node = queue.pop(0)
                    print('======> node: {}'.format(node))
                    component.append(node)
                    for nbr in graph[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
                print('======> visited nodes {}'.format(visited))
                components.append(component)
            else:
                print('-> root {} has been visited'.format(root))
        return components


if __name__ == "__main__":
    graph_adj = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }
    p = ConnectedComponent()

    ccs = list(p.deque(graph_adj))

    print(ccs)

    # print(p.set(graph_adj))