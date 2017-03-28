# nx.yj_layout

This layout applies rule that assign edge length proportional to min(node1.degree, node2.degree) and inverse proportional to abs(node1.degree - node2.degree). This layout tend to strech out a community into several sub-groups with regard to node degrees when two neighboring nodes have high degrees they would be strech out into two groups.