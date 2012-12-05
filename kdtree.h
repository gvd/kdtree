#ifndef KDTREE_H_
#define KDTREE_H_

#include <memory>
#include <limits>
#include <queue>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

template <typename Data, typename Point = boost::geometry::model::d2::point_xy<double>>
class kdtree {
public:
    kdtree() {
        boost::geometry::assert_dimension<Point, 2>();
    }
    virtual ~kdtree() {}
    void add(const Point *point, const Data *data) {
        typename kdnode::ptr node = std::make_shared<kdnode>(point, data);
        m_nodes.push_back(node);
    }
    void build() {
        if (m_nodes.empty()) {
            return;
        }
        m_root = build(m_nodes, 0);
    }
    void clear() {
        m_root.reset();
        m_nodes.clear();
    }
    const Data *nearest_recursive(const Point &query) const {
        if (!m_root) {
            return NULL;
        }
        best_match best(m_root, std::numeric_limits<double>::max());
        nearest(query, m_root, best);
        return best.node->data;
    }
    const Data *nearest_iterative(const Point &query) const {
        if (!m_root) {
            return NULL;
        }
        std::priority_queue<DistanceTuple, std::vector<DistanceTuple>, SmallestOnTop> stack;
        best_match best(m_root, std::numeric_limits<double>::max());
        stack.push(DistanceTuple(0, m_root));
        while (!stack.empty()) {
            const auto current = stack.top();
            if (current.first >= best.distance) {
                return best.node->data;
            }
            stack.pop();
            auto currentNode = current.second;
            double d = boost::geometry::comparable_distance(query, *currentNode->split); // no sqrt
            double dx;
            if (currentNode->axis == 0) {
                dx = boost::geometry::get<0>(query) - boost::geometry::get<0>(*currentNode->split);
            } else {
                dx = boost::geometry::get<1>(query) - boost::geometry::get<1>(*currentNode->split);
            }
            if (d < best.distance) {
                best.node = currentNode;
                best.distance = d;
            }
            node_ptr near = dx <= 0 ? currentNode->left : currentNode->right;
            node_ptr far = dx <= 0 ? currentNode->right : currentNode->left;
            if (far) stack.push(DistanceTuple(dx * dx, far));
            if (near) stack.push(DistanceTuple(0, near));
        }
        return best.node->data;
    }
private:
    struct kdnode {
        typedef std::shared_ptr<kdnode> ptr;
        ptr left;
        ptr right;
        int axis;
        const Point *split;
        const Data *data;
        kdnode(const Point *g, const Data *d) : axis(0), split(g), data(d) {}
    };
    typedef typename kdnode::ptr node_ptr; // get rid of annoying typename
    typedef std::vector<node_ptr> Nodes;
    typedef std::pair<double, node_ptr> DistanceTuple;
    struct SmallestOnTop {
        bool operator()(const DistanceTuple &a, const DistanceTuple &b) const {
            return a.first > b.first;
        }
    };
    Nodes m_nodes;
    node_ptr m_root;

    template<typename NODE_TYPE, std::size_t Dimension>
    struct Sort : std::binary_function<NODE_TYPE, NODE_TYPE, bool> {
        bool operator()(const NODE_TYPE &lhs, const NODE_TYPE &rhs) const {
            return boost::geometry::get<Dimension>(*lhs->split) < boost::geometry::get<Dimension>(*rhs->split);
        }
    };
    struct best_match {
        node_ptr node;
        double distance;
        best_match(const node_ptr &n, double d) : node(n), distance(d) {}
    };

    typename kdnode::ptr build(Nodes &nodes, int depth) {
        if (nodes.empty()) {
            return node_ptr();
        }
        int axis = depth % 2;
        size_t median = nodes.size() / 2;
        if (axis == 0) {
            std::nth_element(nodes.begin(), nodes.begin() + median, nodes.end(), Sort<node_ptr, 0>());
        } else {
            std::nth_element(nodes.begin(), nodes.begin() + median, nodes.end(), Sort<node_ptr, 1>());
        }
        node_ptr node = nodes.at(median);
        node->axis = axis;

        Nodes left(nodes.begin(), nodes.begin() + median);
        Nodes right(nodes.begin() + median + 1, nodes.end());
        node->left = build(left, depth + 1);
        node->right = build(right, depth + 1);

        return node;
    }

    static void nearest(const Point &query, const node_ptr &currentNode, best_match &best) {
        if (!currentNode) {
            return;
        }
        double d = boost::geometry::comparable_distance(query, *currentNode->split); // no sqrt
        double dx;
        if (currentNode->axis == 0) {
            dx = boost::geometry::get<0>(query) - boost::geometry::get<0>(*currentNode->split);
        } else {
            dx = boost::geometry::get<1>(query) - boost::geometry::get<1>(*currentNode->split);
        }
        if (d < best.distance) {
            best.node = currentNode;
            best.distance = d;
        }
        node_ptr near = dx <= 0 ? currentNode->left : currentNode->right;
        node_ptr far = dx <= 0 ? currentNode->right : currentNode->left;
        nearest(query, near, best);
        if ((dx * dx) >= best.distance) {
            return;
        }
        nearest(query, far, best);
    }
};

#endif /* KDTREE_H_ */
