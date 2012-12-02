#include "kdtree.h"

#include <iostream>
#include <vector>

namespace bg = boost::geometry;
namespace bgm = bg::model;

int main(void) {
    typedef bgm::d2::point_xy<double> Point;
    std::vector<Point> points = {{2, 3}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7,2}};

    kdtree<Point, Point> tree;
    for (size_t i = 0; i < points.size(); i++) {
        tree.add(&points[i], &points[i]); // No insert, just adding
    }
    tree.build(); // Bulk build

    Point query(6.5, 1.5);
    const Point *nearest = tree.nearest(query);
    std::cout << "Nearest point to " << bg::dsv(query) << ": " << bg::dsv(*nearest) << std::endl;

    return 0;
}