__author__ = 'Zhan Xianyuan'
# This code is written by Xianyuan Zhan, Purdue University.

# Note:
# Requirement on the node file: the node file must have following format
#                               Node_id, X_Coord, Y_Coord, Invert_Elevation, Area, Outlet_code (1 if the node is a outlet)
# Requirement on the edge file: the first four columns of the edge file must contain follows:
#                               Edge_id, From_Node, To_Node, Length
#                               The edge file should contain at above attribute, these information will be
#                               maintained in the data, however not used in the H-S order computation

import matplotlib.pyplot as plt
import networkx as nx
import math
from copy import deepcopy


# Change the configuration here to run the code
def main():  # For testing purpose
    # Input parameter: parameter to change
    # This is the directory that contains the data
    Fdir = 'data/'
    # This one specify the prefix you want to add to the processed file
    prefix = 'Labeled_'
    # Following are the file names for the node and edge file
    nodeFileName = 'Node_Shintaein_v2.csv'
    edgeFileName = 'Link_Shintaein_v2.csv'

    # Put this one to true if you want to visualize the edge order, for testing only
    # The drawing also shows the direction of the edges
    isDraw = True

    Loopy_HS_order(Fdir, prefix, nodeFileName, edgeFileName, isDraw)


# Major function
def Loopy_HS_order(Fdir, prefix, nodeFileName, edgeFileName, isDraw):
    G = nx.DiGraph()
    # Add nodes to the graph
    root = HS_add_node(Fdir, nodeFileName, G)
    if root == 'NULL':
        print "No root found in the node file"
        quit()

    # print root
    # Add edges to the graph
    HS_add_edge(Fdir, edgeFileName, G)

    # We have to work on a network that every node are connected to the root, if not we need to prune the network
    HS_get_connected_component(G, root)

    print 'Number of nodes loaded: ' + str(G.number_of_nodes())
    print 'Number of edges loaded: ' + str(G.number_of_edges())

    # Reverse the graph in order to compute post dominators
    revG = G.reverse(True)
    # Here we compute the post dominator tree of the graph
    doms = nx.immediate_dominators(revG, root)
    # remove teh reverse graph to save memory
    revG.clear()

    in_nodes = HS_get_leaves(G)
    print 'Number of initial in-nodes = ' + str(len(in_nodes))

    in_nodes = HS_simplify_from_top(G, in_nodes)

    # This is main iteration
    iteration = 1
    while 1:
    # for i in range(3):
        print "\nIteration: " + str(iteration)
        iteration += 1

        end_nodes = HS_resolve_super_closure(G, in_nodes, doms)

        if len(end_nodes) == 0:
            break
        else:
            in_nodes = end_nodes

    # Work done! Write results to files
    HS_output_graph(G, Fdir, prefix, nodeFileName, edgeFileName)
    if isDraw:
        HS_draw_graph(G)


def HS_add_node(Fdir, nodeFileName, G):
    # First we need to parse the node file
    nodeFile = open(Fdir + nodeFileName, 'r')

    # Escape the header
    line = nodeFile.readline()
    root = 'NULL'

    while 1:
        line = nodeFile.readline()
        if len(line) == 0:
            break

        node_id = line.split(',')[0]
        x_coord = float(line.split(',')[1])
        y_coord = float(line.split(',')[2])
        invert_elevation = float(line.split(',')[3])
        node_area = float(line.split(',')[4])
        flag = int(line.split(',')[5].rstrip('\n'))

        if flag == 1:
            # return the node key of the outlet
            root = node_id
        # Add the node
        # G.add_node(node_id, pos=(x_coord, y_coord))
        G.add_node(node_id, pos=(x_coord, y_coord), invElevation=invert_elevation, area=node_area, cumArea=0, used=0)

    nodeFile.close()
    return root


def HS_add_edge(Fdir, edgeFileName, G):
    # First parse the edge file
    edgeFile = open(Fdir + edgeFileName, 'r')

    # Escape the header
    line = edgeFile.readline()
    header = line.rstrip('\n').split(',')
    Nattr = len(line.split(','))
    if Nattr < 4:
        print "Problematic edge file!"
        quit()

    while 1:
        line = edgeFile.readline()
        if len(line) == 0:
            break

        line = line.rstrip('\n')
        edge_id = line.split(',')[0]
        from_node = line.split(',')[1]
        to_node = line.split(',')[2]
        node_length = float(line.split(',')[3])

        # Rest of the attribute will be put into a separate dictionary
        edge_otherAttr = {}
        edge_otherAttr['edge_id'] = edge_id
        if Nattr > 4:
            for i in range(Nattr - 4):
                edge_otherAttr[header[4 + i]] = line.split(',')[4 + i]

        # Add the edges
        G.add_edge(from_node, to_node, length=node_length, used=0, order=0, otherAttr=edge_otherAttr)

    edgeFile.close()


# Following functions are the algorithms part:
# Get a sub-graph in which every nodes are connected to the root
def HS_get_connected_component(G, root):
    # These are all the nodes that have path to root
    reachable_nodes = nx.ancestors(G, root)
    # print reachable_nodes
    # We iterate through all the nodes and remove those that are not reachable
    for node_key in G.nodes():
        if node_key not in reachable_nodes:
            if node_key != root:
                G.remove_node(node_key)


# Identify all the leaves and initialize the peripheral nodes
def HS_get_leaves(G):
    leaves = []
    in_nodes = []
    for node_key in G.nodes():
        if G.in_degree(node_key) == 0:
            if G.out_degree(node_key) > 0:
                leaves.append(node_key)
                # Marke the leaves as used and initialize the cumulative area
                G.node[node_key]['used'] = 1
                G.node[node_key]['cumArea'] = G.node[node_key]['area']

                # Next we need to mark the initial first order stream edges
                successors = G.successors(node_key)
                # print successors
                for child in successors:
                    G[node_key][child]['used'] = 1
                    G[node_key][child]['order'] = 1
                    if child not in in_nodes:
                        in_nodes.append(child)

    return in_nodes


# This implements the Horton-Strauler order on tree like structures
# We define a recursive procedure to handle the task
# The function always returns a set of unconcretized node with concretized link end at it
def HS_simplify_from_top(G, in_nodes):
    # For each node, we do one step further, and stop only if for all nodes, there is no improvement
    # Number of nodes improved
    N_improve = 0
    for node in in_nodes:
        # We first get its out-going neighbors
        successors = G.successors(node)
        predecessors = G.predecessors(node)

        if len(successors) == 1:
            child = successors[0]
            # We have encountered a tree like structure, we can potentially improve it
            # We first need to check if all its predecessors are concretized, if so we can directly improve,
            # otherwise, we need to wait until all the other ancestors are concretized
            N_concretized = 0
            for parent in predecessors:
                if G[parent][node]['used'] == 1:
                    N_concretized += 1
            if N_concretized == len(predecessors):
                upstream_order = []
                upstream_area = []
                for parent in predecessors:
                    upstream_order.append(G[parent][node]['order'])
                    upstream_area.append(G.node[parent]['cumArea'])

                # Determine the Horton-Strauler stream order
                G[node][child]['order'] = HS_order_rule(upstream_order)
                G[node][child]['used'] = 1
                N_improve += 1

                # replace the original in-node with the child
                for i in range(len(in_nodes)):
                    if in_nodes[i] == node:
                        G.node[node]['used'] = 1
                        G.node[node]['cumArea'] = sum(upstream_area)
                        in_nodes[i] = child
                        break

        else:
            if len(successors) == 0:
                # This stream died, although unlikely, just process and remove the node from in-node list
                G.node[node]['used'] = 1
                upstream_area = []
                for parent in predecessors:
                    upstream_area.append(G.node[parent]['cumArea'])
                G.node[node]['cumArea'] = sum(upstream_area)
                # Remove the node
                in_nodes.remove(node)
                N_improve += 1

    # We need to filter the duplicates
    in_nodes = list(set(in_nodes))
    # OK, here is the recursion part, stop only when no improvement can be make
    if N_improve > 0:
        in_nodes = HS_simplify_from_top(G, in_nodes)
        return in_nodes
    else:
        return in_nodes


# Original Horton-Strauler stream order rule when performing simplify-from-top
def HS_order_rule(upstream_order):
    max_order = max(upstream_order)
    # count number of max_order stream
    N_max = 0
    for order in upstream_order:
        if order == max_order:
            N_max += 1

    return HS_merge_order_rule(N_max, max_order)


# This parts implement HS stream order merging rule in the undirected loopy sub-structure
def HS_merge_order_rule(N_max, max_order):
    if N_max > 1:
        return max_order + 1
    else:
        return max_order

# # This one is the old function: Obsolete
# # This one handles the super-closure, which performs simplify from top, and forms sub-closures, and finally returns IPDs as new in-nodes
# def HS_resolve_super_closure(G, in_nodes, doms):
#     # Initial simplify_from_top
#     in_nodes = HS_simplify_from_top(G, in_nodes)
#
#     print "Number of in_nodes = " + str(len(in_nodes))
#     print in_nodes
#     # From the in_nodes, we form a set of sub-closures
#     sub_closure = {}
#     IPDdict = {}
#     for node_key in in_nodes:
#         IPD = doms[node_key]
#         # We might encounter cases that the in-nodes itself serve as IPD, under such cases, we need to merge it into downstream dominators
#         while 1:
#             if IPD in in_nodes:
#                 IPD = doms[IPD]
#             else:
#                 break
#         if IPD not in sub_closure.keys():
#             sub_closure[IPD] = [node_key]
#         else:
#             sub_closure[IPD].append(node_key)
#         IPDdict[node_key] = IPD
#
#     # We need to continue merge some of the sub-closure, that the in-nodes have in-nodes ancestors
#     # Under such cases, we merge the sub-closure to its in-nodes ancestor's sub-closure
#     flag = True
#     while flag:
#         flag = False
#         for node_key in in_nodes:
#             predecessors = G.predecessors(node_key)
#             local_in_node = set(predecessors).intersection(sub_closure[IPDdict[node_key]])
#             all_in_node = set(predecessors).intersection(in_nodes)
#             out_in_node = list(all_in_node - local_in_node)
#             if len(out_in_node) >0:
#                 # We are in trouble, we need to merge the sub-closure to its parents'
#                 # We just merge to the first parent, the rest will be handled by the loop
#                 # directly merge
#                 sub_closure[IPDdict[out_in_node[0]]] = sub_closure[IPDdict[out_in_node[0]]] + sub_closure[IPDdict[node_key]]
#                 old_IPD = IPDdict[node_key]
#                 for node in sub_closure[old_IPD]:
#                     IPDdict[node] = IPDdict[out_in_node[0]]
#                 del sub_closure[old_IPD]
#                 flag = True
#
#     print "Number of IPDs = " + str(len(sub_closure.keys()))
#     N_improve = 0
#     IPD_list = sub_closure.keys()
#     for IPD in IPD_list:
#         print "Sub-closure: " + IPD
#         print sub_closure[IPD]
#
#         solvable = HS_solvability_check(G, sub_closure[IPD], IPD)
#         if solvable:
#             end_nodes = HS_resolve_sub_closure(G, sub_closure[IPD], IPD)
#             # replace the in_nodes with the IPD in the super_closure
#             for i in range(len(in_nodes)):
#                 if in_nodes[i] in sub_closure[IPD]:
#                     in_nodes[i] = end_nodes[0]
#
#             if len(end_nodes)>1:
#                 for i in range(len(end_nodes)-1):
#                     in_nodes.append(end_nodes[i+1])
#             # Remove this sub-closure entry
#             del sub_closure[IPD]
#             N_improve += 1
#
#     # # If we got stuck, then we need to further analyze the dependency among sub-closures
#     # IPD_dependent = {}
#     # for IPD in sub_closure.keys():
#     #     IPD_dependent[IPD] = []
#
#     # if N_improve == 0:
#     #     for node_key in in_nodes:
#     #         ancestor_nodes = nx.ancestors(G,node_key)
#     #         local_in_node = ancestor_nodes.intersection(sub_closure[IPDdict[node_key]])
#     #         all_in_node = ancestor_nodes.intersection(in_nodes)
#     #         out_in_node = list(all_in_node - local_in_node)
#     #         if len(out_in_node) >0:
#     #             for node in out_in_node:
#     #                 if IPDdict[node] not in IPD_dependent[IPDdict[node_key]]:
#     #                     IPD_dependent[IPDdict[node_key]].append(IPDdict[node])
#     #     print IPD_dependent
#     #
#     #     # We merge mutually dependent sub-closures, here we just do 1 if necessary:
#     #     # if a depend on b, b depend on a, then new IPD = IPD(IPD(a), IPD(b))
#     #     IPD_merge_list = []
#     #     for IPD in sub_closure.keys():
#     #         if IPD not in IPD_merge_list:
#     #             for depIPD in IPD_dependent[IPD]:
#     #                 if IPD in IPD_dependent[depIPD]:
#     #                     IPD_merge_list.append(IPD)
#     #                     IPD_merge_list.append(depIPD)
#
#     in_nodes = list(set(in_nodes))
#
#     # return ['Out']  # for testing only
#     return in_nodes

# This one handles the super-closure, which performs simplify from top, and forms sub-closures, and finally returns IPDs as new in-nodes
def HS_resolve_super_closure(G, in_nodes, doms):
    # Initial simplify_from_top
    in_nodes = HS_simplify_from_top(G, in_nodes)

    print "Number of in_nodes = " + str(len(in_nodes))
    print in_nodes
    # From the in_nodes, we form a set of sub-closures
    sub_closure = {}
    IPDdict = {}
    for node_key in in_nodes:
        IPD = doms[node_key]
        # We might encounter cases that the in-nodes itself serve as IPD, under such cases, we need to merge it into downstream dominators
        while 1:
            if IPD in in_nodes:
                IPD = doms[IPD]
            else:
                break
        if IPD not in sub_closure.keys():
            sub_closure[IPD] = [node_key]
        else:
            sub_closure[IPD].append(node_key)
        IPDdict[node_key] = IPD

    print "Initial Number of IPDs = " + str(len(sub_closure.keys()))
    print sub_closure.keys()
    # We need to analyze the dependency between IPDs, so as to merge the mutually dependent IPDs
    IPD_dependent = {}
    for IPD in sub_closure.keys():
        IPD_dependent[IPD] = []

    for node_key in in_nodes:
        ancestor_nodes = nx.ancestors(G,node_key)
        local_in_node = ancestor_nodes.intersection(sub_closure[IPDdict[node_key]])
        all_in_node = ancestor_nodes.intersection(in_nodes)
        out_in_node = list(all_in_node - local_in_node)
        if len(out_in_node) >0:
            for node in out_in_node:
                if IPDdict[node] not in IPD_dependent[IPDdict[node_key]]:
                    IPD_dependent[IPDdict[node_key]].append(IPDdict[node])
    print IPD_dependent

    IPD_merge_list = []
    for IPD in sub_closure.keys():
        # The first case, if the sub-closure does not depend on anything, we can solve directly
        if len(IPD_dependent[IPD]) == 0:
            IPD_merge_list.append([IPD])
        else:
            # The sub-closure depend on other sub-closures
            cand_set = [IPD]
            for depIPD in IPD_dependent[IPD]:
                if IPD in IPD_dependent[depIPD]:
                    # if it mutually dependent with another IPD
                    cand_set.append(depIPD)
            # We now have a candidate merge set, we need to put it into a proper entry in IPD_merge_list
            flag = True
            for i in range(len(IPD_merge_list)):
                if len(set(IPD_merge_list[i]).intersection(cand_set)) > 0:
                    flag = False
                    # We've found an entry to replace
                    for node in cand_set:
                        if node not in IPD_merge_list[i]:
                            IPD_merge_list[i].append(node)
                    break
            if flag:
                 IPD_merge_list.append(cand_set)

    print "Mergable IPDs:"
    print IPD_merge_list

    # Now we need to merge these mergable IPDs
    for merge_list in IPD_merge_list:
        if len(merge_list) > 1:
            # We need to do some merging here
            # We use the IPDs in the merge_list as in-nodes to create a newIPD
            temp_in_nodes = deepcopy(merge_list)
            # We first attampt to reduce the set by only keeping the dominating IPDs
            # We need to assume acyclic graph here, otherwise we are in trouble
            for nodeE in merge_list:
                for nodeS in merge_list:
                    if nodeS != nodeE:
                        if nx.has_path(G, nodeS, nodeE):
                            # remove nodeS
                            if nodeS in temp_in_nodes:
                                temp_in_nodes.remove(nodeS)
            # print temp_in_nodes
            if len(temp_in_nodes) == 1:
                print merge_list
                # We find a singleton, we merge everything together
                cand_closure = sub_closure[temp_in_nodes[0]]
                for old_IPD in merge_list:
                    if old_IPD != temp_in_nodes[0]:
                        cand_closure = cand_closure + sub_closure[old_IPD]
                        del sub_closure[old_IPD]
                sub_closure[temp_in_nodes[0]] = cand_closure

    # print "Number of IPDs = " + str(len(sub_closure.keys()))
    IPD_list = sub_closure.keys()
    for IPD in IPD_list:
        print "Sub-closure: " + IPD
        print sub_closure[IPD]

        solvable = HS_solvability_check(G, sub_closure[IPD], IPD)
        if solvable:
            end_nodes = HS_resolve_sub_closure(G, sub_closure[IPD], IPD)
            # replace the in_nodes with the IPD in the super_closure
            for i in range(len(in_nodes)):
                if in_nodes[i] in sub_closure[IPD]:
                    in_nodes[i] = end_nodes[0]

            if len(end_nodes)>1:
                for i in range(len(end_nodes)-1):
                    in_nodes.append(end_nodes[i+1])
            # Remove this sub-closure entry
            del sub_closure[IPD]

    in_nodes = list(set(in_nodes))

    return in_nodes

# Function checks if a sub-closure is solvable
def HS_solvability_check(G, in_nodes, IPD):
    solvable = True
    sub_nodes = []
    for node in in_nodes:
        sub_nodes = HS_get_decedents_before(G, sub_nodes, node, IPD)

    # For each in-nodes, we check if all its external in-links are concretized
    for node in sub_nodes:
        predecessors = G.predecessors(node)
        if len(predecessors) > 0:
            # We should not consider a in-node if it is not the out-most one
            ext_predecessors = set(predecessors) - set(sub_nodes)
            for parent in ext_predecessors:
                if G[parent][node]['used'] == 0:
                    # we allow unconcretized in_nodes, since they are not the out-most in-nodes
                    print "IPD: parent -- node: " + IPD + ': ' + parent + ' - ' + node
                    solvable = False
                    break

    return  solvable



def HS_resolve_sub_closure(G, in_nodes, IPD):
    # Since the out-most in-nodes' parents are all concretized, thus we can safely extract the sub-graph between in_nodes and IPD
    sub_nodes = []
    for node in in_nodes:
        sub_nodes = HS_get_decedents_before(G, sub_nodes, node, IPD)

    print "Solvable Sub-closure: " + IPD
    print in_nodes
    print sub_nodes

    # Extract the sub-graph
    subGraph = G.subgraph(sub_nodes)

    # Create the initial order set
    orderSet = {}

    # initialize the order at the edge at IPD
    # We first check if N_in^max >1
    O_in_max = 0
    O_in_min = 999
    N_in_max = 0
    in_node_info = {}
    in_node_upstream_area = []
    for node in in_nodes:
        # in-nodes at least have 1 concrete parent
        in_node_info[node] = HS_get_anc_info(G, node) # ancestors of in-nodes not in subGraph
        if O_in_max < in_node_info[node]["max_order"]:
            O_in_max = in_node_info[node]["max_order"]
            N_in_max = 1
        else:
            if O_in_max == in_node_info[node]["max_order"]:
                N_in_max += 1
        if O_in_min > in_node_info[node]["min_order"]:
            O_in_min = in_node_info[node]["min_order"]
        in_node_upstream_area = in_node_upstream_area + in_node_info[node]["upstream_area"]

    end_nodes = G.successors(IPD) # decendents of end-nodes out of IPD
    end_edges = [] # These edges are also not in subGraph
    for enode in end_nodes:
        end_edges.append((IPD, enode))

    if N_in_max > 1:
        # For all edges out of IPD, we assign order starting from O_in_max + 1
        for edge in end_edges:
            lower = O_in_max + 1
            upper = O_in_max + int(math.floor(math.log(len(in_nodes), 2)))

            if upper > lower:
                orderSet[edge] = range(lower, upper + 1)
            else:
                orderSet[edge] = [min(lower, upper)]
                # order set is a singleton, we can directly concretize it
                HS_concretize(G, edge, orderSet[edge][0], in_node_upstream_area)
    else:
        # Otherwise, we start the order from O_in_max
        for edge in end_edges:
            lower = O_in_max
            upper = O_in_max + int(math.floor(math.log(len(in_nodes), 2)))

            if upper > lower:
                orderSet[edge] = range(lower, upper + 1)
            else:
                orderSet[edge] = [min(lower, upper)]
                # order set is a singleton, we can directly concretize it
                HS_concretize(G, edge, orderSet[edge][0], in_node_upstream_area)

    # Next, we initialize the OrderSet of the rest of the edges
    for edge in subGraph.edges():
        lower = O_in_min
        upper = O_in_max + int(math.floor(math.log(len(in_nodes), 2)))

        if upper > lower:
            orderSet[edge] = range(lower, upper + 1)
        else:
            orderSet[edge] = [min(lower, upper)]
            # order set is a singleton, we can directly concretize it
            HS_concretize(G, edge, orderSet[edge][0], [0])

    # orderSet construction finished, we don't need the subGraph anymore, we can remove it
    subGraph.clear()

    # # debugging:
    # for edge in orderSet.keys():
    #     print edge
    #     print orderSet[edge]

    # Number of unconcreteized edges
    N_unconcrete = len(orderSet.keys())

    iteration = 0
    while N_unconcrete > 0:
        iteration += 1
        N_improve = 0
        # First handle the diverging case
        for edge in orderSet.keys():
            if G[edge[0]][edge[1]]['used'] == 0:
                # Let's get the ancestor's info to see if it is all concretized
                anc_info = HS_get_anc_info(G, edge[0])
                if anc_info["N_concrete"] == anc_info["N_parents"]:
                    # check if it has diverges
                    successors = G.successors(edge[0])
                    if len(successors) > 0:
                        for node in successors:
                            orderSet[(edge[0], node)] = [anc_info["max_order"]]
                            print "Diverging handled: edge: " + edge[0] + ' - ' + edge[1] + " max_order=" + str(anc_info["max_order"])
                            HS_concretize(G, (edge[0], node), anc_info["max_order"], [0])
                            N_improve += 1

        # We count the exact number of N_unconcrete, this step perform together with prune order set
        cur_N_unconcrete = 0
        for edge in orderSet.keys():
            if G[edge[0]][edge[1]]['used'] == 0:
                flag = HS_prune_order_set(G, edge, orderSet)
                if flag < 2:
                    cur_N_unconcrete += 1
                if flag > 0:
                    N_improve += 1
                    print "Edge pruned: edge: " + edge[0] + ' - ' + edge[1]

        # This parts handles the merging case
        if N_improve == 0 and cur_N_unconcrete > 0:
            for edge in orderSet.keys():
                if G[edge[0]][edge[1]]['used'] == 0:
                    # We need to check if it is a merging with all concretized ancestors
                    anc_info = HS_get_anc_info(G, edge[0])
                    if G.out_degree(edge[0]) == 1:
                        if anc_info["N_concrete"] == anc_info["N_parents"] and anc_info["N_parents"] > 1:
                            # This is a valid merging, we can merge it
                            orderSet[edge] = [HS_merge_order_rule(anc_info["N_max_order"], anc_info["max_order"])]
                            HS_concretize(G, edge, orderSet[edge][0], [0])
                            N_improve += 1
                            print "Improved merging: edge: " + edge[0] + ' - ' + edge[1]
                            print G.out_degree(edge[0])
                            # We just do one merging at a time
                            break

        # Update the number of unconcretized edges
        if cur_N_unconcrete < N_unconcrete:
            N_unconcrete = cur_N_unconcrete


        if N_unconcrete > 0 and N_improve == 0:
            print "I can't solve this sub-closure any further"
            # We exit early, so we need to revise the end_nodes to include all nodes that are not concretized
            end_nodes = []
            for edge in orderSet.keys():
                if G[edge[0]][edge[1]]['used'] == 0:
                    anc_info = HS_get_anc_info(G, edge[0])
                    if anc_info["N_concrete"] > 0:
                        if edge[0] not in end_nodes:
                            end_nodes.append(edge[0])
            break


    return end_nodes

# concretize a edge and its from node
# Set upstream_area to [0] if do not need to compute
def HS_concretize(G, edge, edgeOrder, upstream_area):
    G[edge[0]][edge[1]]['order'] = edgeOrder
    G[edge[0]][edge[1]]['used'] = 1
    G.node[edge[0]]['used'] = 1
    G.node[edge[0]]['cumArea'] = sum(upstream_area)
    print "Concretized: " + edge[0] + " - " + edge[1] + " Order=" + str(edgeOrder)

# return the information related to the stream orders of the ancestor
def HS_get_anc_info(G, node):
    anc_info = {}
    predecessors = G.predecessors(node)
    anc_info["N_parents"] = len(predecessors)
    anc_info["N_concrete"] = 0
    anc_info["N_max_order"] = 0
    anc_info["min_order"] = 999
    anc_info["max_order"] = 0
    anc_info["upstream_area"] = []

    for parent in predecessors:
        if G[parent][node]['order'] >0:
            # we find a concretized parent
            anc_info['N_concrete'] += 1
            if G[parent][node]['order'] > anc_info["max_order"]:
                anc_info["N_max_order"] = 1
                anc_info["max_order"] = G[parent][node]['order']
            else:
                if G[parent][node]['order'] == anc_info["max_order"]:
                    anc_info["N_max_order"] += 1

            if G[parent][node]['order'] < anc_info["min_order"]:
                anc_info["min_order"] = G[parent][node]['order']

        # This one gets the list of upstream_areas for computing downstream cumArea
        if G.node[parent]['cumArea'] > 0:
            anc_info["upstream_area"].append(G.node[parent]['cumArea'])

    return anc_info

# return the information related to the stream orders of the decendents
def HS_get_suc_info(G, node):
    suc_info = {}
    successors = G.successors(node)
    suc_info["N_children"] = len(successors)
    suc_info["N_concrete"] = 0
    suc_info["min_order"] = 999
    suc_info["max_order"] = 0

    for child in successors:
        if G[node][child]['order'] >0:
            # we find a concretized parent
            suc_info['N_concrete'] += 1
            if G[node][child]['order'] > suc_info["max_order"]:
                suc_info["max_order"] = G[node][child]['order']

            if G[node][child]['order'] < suc_info["min_order"]:
                suc_info["min_order"] = G[node][child]['order']

    return suc_info


# This part uses rules to prune order set.
# Return following information: 0- no improvement; 1- improved, unconcretized; 2- improved and concretized
def HS_prune_order_set(G, edge, orderSet):
    candOrder = orderSet[edge]
    nOrder = len(candOrder)
    flag = 0

    if nOrder == 1:
        HS_concretize(G, edge, candOrder[0], [0])
        return 2

    lower = candOrder[0]
    upper = candOrder[nOrder - 1]

    # get ancestor and successor's info
    anc_info = HS_get_anc_info(G, edge[0])
    suc_info = HS_get_suc_info(G, edge[1])

    # Rule 1: update lower
    if anc_info["N_concrete"] > 0 and lower < anc_info["max_order"]:
    # Prune out all orders that below O_anc_max
        for order in candOrder:
            if order < anc_info["max_order"]:
                candOrder.remove(order)
                nOrder -= 1
        lower = candOrder[0]
        # flag = 1

    # Rule 2: update the upper
    new_upper = upper
    if anc_info["N_concrete"] > 0 and anc_info["N_max_order"] > 1:
        if upper > anc_info["max_order"] + 1:
            new_upper = anc_info["max_order"] + 1

    if suc_info["N_concrete"] > 0 and new_upper > suc_info["min_order"]:
        new_upper = suc_info["min_order"]

    if upper > new_upper:
        for order in candOrder:
            if order > new_upper:
                candOrder.remove(order)
                nOrder -= 1
        # flag = 1

    upper = candOrder[nOrder-1]

    # Rule 3: additional rule: no-merge
    # If it is a merge
    if G.in_degree(edge[1]) > 1 and nOrder > 0:
        # get successor's ancestor info
        dec_anc_info = HS_get_anc_info(G, edge[1])
        if anc_info["N_concrete"] > 0:
            if lower < suc_info["min_order"] and suc_info["min_order"] == dec_anc_info["max_order"]:
                upper = suc_info["min_order"] - 1
                for order in candOrder:
                    if order > upper:
                        candOrder.remove(order)
                        nOrder -= 1
                # flag = 1


    if candOrder < orderSet[edge]:
        orderSet[edge] = candOrder
        flag = 1
    else:
        flag = 0

    # We concretize the edge if it is a singleton
    if nOrder == 1:
        HS_concretize(G, edge, candOrder[0], [0])
        flag = 2

    return flag


# Get all decendents (including the node itself) from node to end node
def HS_get_decedents_before(G, node_list, from_node, end_node):
    if from_node not in node_list:
        node_list.append(from_node)
    successors = G.successors(from_node)
    if len(successors) == 1 and successors[0] == end_node:
        if end_node not in node_list:
            node_list.append(end_node)
        return node_list

    for node in successors:
        node_list = HS_get_decedents_before(G, node_list, node, end_node)

    return node_list


# A drawing function, used for debugging
def HS_draw_graph(G):
    pos = nx.get_node_attributes(G, 'pos')
    order_label = nx.get_edge_attributes(G, 'order')
    plt.figure(figsize=(10, 10))
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, width=0.8)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=order_label, label_pos=0.5, font_size=10, alpha=1, rotate=False)
    nx.draw_networkx_labels(G, pos, font_size=15)

    plt.savefig('test.png')
    plt.show()


# Output data
def HS_output_graph(G, Fdir, prefix, nodeFileName, edgeFileName):
    print "Writing processed graph to file..."
    nodeFile = open(Fdir + prefix + nodeFileName, 'w')
    edgeFile = open(Fdir + prefix + edgeFileName, 'w')

    # Writing to the node file
    outputline = "Node_id,X-coord,Y-Coord,Invert_Elevation,Area,Cumulative_Area\n"
    nodeFile.write(outputline)
    for node in G.nodes():
        outputline = node + ',' + str(G.node[node]['pos'][0]) + ',' + str(G.node[node]['pos'][1]) + \
                     ',' + str(G.node[node]['invElevation']) + ',' + str(G.node[node]['area']) + ','
        cArea = G.node[node]['cumArea']
        if G.node[node]['cumArea'] == 0:
            cArea = -1
        outputline = outputline + str(cArea) + '\n'
        nodeFile.write(outputline)

    nodeFile.close()

    # Writing to edge file
    outputline = "Link_ID,From_Node,To_Node,Length,Slope,Stream_Order"
    for header_key in G[G.edges()[0][0]][G.edges()[0][1]]['otherAttr'].keys():
        if header_key != "edge_id":
            outputline = outputline + ',' + header_key
    outputline = outputline + '\n'
    edgeFile.write(outputline)
    for edge in G.edges():
        elevation1 = G.node[edge[0]]['invElevation']
        elevation2 = G.node[edge[1]]['invElevation']
        slope = abs(elevation1 - elevation2) / G[edge[0]][edge[1]]['length']
        outputline = G[edge[0]][edge[1]]['otherAttr']['edge_id'] + ',' + edge[0] + ',' + edge[1] \
                     + ',' + str(G[edge[0]][edge[1]]['length']) + ',' + str(slope) + ',' \
                     + str(G[edge[0]][edge[1]]['order'])
        for header_key in G[edge[0]][edge[1]]['otherAttr'].keys():
            if header_key != "edge_id":
                outputline = outputline + ',' + G[edge[0]][edge[1]]['otherAttr'][header_key]
        outputline = outputline + '\n'
        edgeFile.write(outputline)

    edgeFile.close()

if __name__ == "__main__":
    main()
