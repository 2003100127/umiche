__version__ = "v1.0"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"


from typing import List, Dict

from collections import deque
from umiche.util.Console import Console


class CC:

    def __init__(
            self,
            verbose=False,
    ):
        self.console = Console()
        self.console.verbose = verbose


    def deque(
            self,
            graph : Dict,
    ):
        """

        Parameters
        ----------
        graph
            1d dict of an adjacency list

        Returns
        -------
        2d list, each 1d list representing a component

        """
        visited = set()
        for root, nbrs in graph.items():
            if root not in visited:
                visited.add(root)
                component = []
                queue = deque([root])
                self.console.print('======> root {} has not been visited'.format(root))
                self.console.print('======> a queue built by root {} is {}'.format(root, queue))
                while queue:
                    self.console.print('=========> a queue built by each root node {}'.format(queue))
                    node = queue.popleft()
                    self.console.print('=========> node: {}'.format(node))
                    component.append(node)
                    for nbr in graph[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
                self.console.print('=========> visited nodes {}'.format(visited))
                yield component
            else:
                self.console.print('=========> root {} has been visited'.format(root))
                continue

    def set(
            self,
            graph : Dict,
    ) -> List:
        """
        Examples
        --------
        >>>seen = set()
        >>>def component(node):
        >>>    nodes = set([node])
        >>>    while nodes:
        >>>        node = nodes.pop()
        >>>        seen.add(node)
        >>>        nodes |= neighbors[node] - seen
        >>>        yield node
        >>>    for node in neighbors:
        >>>        if node not in seen:
        >>>            yield component(node)

        Parameters
        ----------
        graph
            1d dict of an adjacency list

        Returns
        -------
        2d list, each 1d list representing a component

        """
        visited = set()
        components = []
        for root, nbrs in graph.items():
            if root not in visited:
                visited.add(root)
                component = []
                queue = [root]
                self.console.print('======> root {} has not been visited'.format(root))
                self.console.print('======> a queue built by root {} is {}'.format(root, queue))
                while queue:
                    self.console.print('======> a queue built by each root node {}'.format(queue))
                    node = queue.pop(0)
                    self.console.print('======> node: {}'.format(node))
                    component.append(node)
                    for nbr in graph[node]:
                        if nbr not in visited:
                            visited.add(nbr)
                            queue.append(nbr)
                self.console.print('======> visited nodes {}'.format(visited))
                components.append(component)
            else:
                self.console.print('=========>root {} has been visited'.format(root))
        return components


if __name__ == "__main__":
    graph_adj_mclumi = {
        'A': ['B', 'C', 'D'],
        'B': ['A', 'C'],
        'C': ['A', 'B'],
        'D': ['A', 'E', 'F'],
        'E': ['D', 'G'],
        'F': ['D', 'G'],
        'G': ['E', 'F'],
    }
    p = CC()

    print(list(p.deque(graph=graph_adj_mclumi)))

    print(p.set(graph=graph_adj_mclumi))