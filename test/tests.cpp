#include <gtest/gtest.h>

#include <vector>
#include <iostream>
#include <chrono>

#include <boost/random.hpp> // mt19937
#include <boost/nondet_random.hpp> // random_device
#include <boost/random/uniform_real_distribution.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <kdtree.h>

namespace bg = boost::geometry;
namespace bgm = bg::model;
namespace bacc = boost::accumulators;

typedef bgm::d2::point_xy<double> Point;

typedef bacc::accumulator_set < size_t, bacc::stats <bacc::tag::mean, bacc::tag::min,
    bacc::tag::max, bacc::tag::variance, bacc::tag::median, bacc::tag::sum, bacc::tag::count>> TimeAccumulator;

std::ostream &operator<<(std::ostream &out, const TimeAccumulator &acc) {
    out << "count: " << bacc::count(acc) << "\n";
    out << "min: " << bacc::min(acc) << "\n";
    out << "max: " << bacc::max(acc)  << "\n";
    out << "mean: " << bacc::mean(acc) << "\n";
    out << "variance: " << bacc::variance(acc) << "\n";
    out << "median: " << bacc::median(acc) << "\n";
    out << "sum: " << bacc::sum(acc);
    return out;
}

template <typename Point>
void randomPoints(size_t nr, std::vector<Point> &points) {
    points.resize(nr);
    boost::mt19937 gen(time(0));
    boost::random::uniform_real_distribution<> dis_x(-10 , 10);
    boost::random::uniform_real_distribution<> dis_y(-10 , 10);
    for (size_t i = 0; i < nr; i++) {
        Point p = {dis_x(gen), dis_y(gen)};
        points.push_back(p);
    }
}

class KdTreeTest : public ::testing::Test {
public:
    KdTreeTest() : m_gen(time(0)), m_dis_x(-10, 10),
        m_dis_y(-10, 10), m_query_count(1000), m_point_count(100000) {
    }
protected:
     virtual void SetUp() {
        randomPoints(m_point_count, m_points);
        for (size_t i = 0; i < m_points.size(); i++) {
            m_tree.add(&m_points[i], &m_points[i]); // No insert, just adding
        }
        m_tree.build(); // Bulk build
    }
    virtual void TearDown() {
        m_tree.clear();
        m_points.clear();
    }

    Point RandomPoint() {
        return {m_dis_x(m_gen), m_dis_y(m_gen)};
    }
    std::vector<Point> m_points;
    kdtree<Point> m_tree;
    boost::mt19937 m_gen;
    boost::random::uniform_real_distribution<> m_dis_x;
    boost::random::uniform_real_distribution<> m_dis_y;
    size_t m_query_count;
    size_t m_point_count;
};

TEST_F(KdTreeTest, identical_iterative_and_recursive_results) {
    for (size_t i = 0; i < m_query_count; i++) {
        const Point query = RandomPoint();
        auto recursive_result = m_tree.nearest_recursive(query);
        auto iterative_result = m_tree.nearest_iterative(query);
        bool identical = bg::equals(*recursive_result, *iterative_result);
        EXPECT_TRUE(recursive_result != NULL);
        EXPECT_TRUE(identical) << bg::dsv(*recursive_result) << " != "<< bg::dsv(*iterative_result);
    }
}

TEST_F(KdTreeTest, recursive_performance) {
    TimeAccumulator time_acc;
    for (size_t i = 0; i < m_query_count; i++) {
        const Point query = RandomPoint();
        auto startTime = std::chrono::high_resolution_clock::now();
        const Point *nearest = m_tree.nearest_recursive(query);
        auto endTime = std::chrono::high_resolution_clock::now();
        size_t t = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        time_acc(t);
    }
    std::cout << "Recursive performance (usec):\n" << time_acc << std::endl;
}

TEST_F(KdTreeTest, iterative_performance) {
    TimeAccumulator time_acc;
    for (size_t i = 0; i < m_query_count; i++) {
        const Point query = RandomPoint();
        auto startTime = std::chrono::high_resolution_clock::now();
        const Point *nearest = m_tree.nearest_iterative(query);
        auto endTime = std::chrono::high_resolution_clock::now();
        size_t t = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        time_acc(t);
    }
    std::cout << "Iterative performance (usec):\n" << time_acc << std::endl;
}

TEST(WikipediaExample, test) {
    std::vector<Point> points = {{2, 3}, {5, 4}, {9, 6}, {4, 7}, {8, 1}, {7,2}};
    kdtree<Point> tree;
    for (size_t i = 0; i < points.size(); i++) {
        tree.add(&points[i], &points[i]);
    }
    tree.build();
    const Point query(5, 6);
    const Point *nearest = tree.nearest_recursive(query);
    const Point expected_result(4, 7);
    bool identical = bg::equals(*nearest, expected_result);
    EXPECT_TRUE(nearest != NULL);
    EXPECT_TRUE(identical) << bg::dsv(*nearest) << " != " << bg::dsv(expected_result);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
