#pragma once

#include <vector>

void test_PowerSpec();

void test_vec();

void test_profile_solver();

double test_func(double x, void* params);

void test_boostmath();

void test_gsl_integration();

namespace PowerTest {

double mean(const std::vector<double>& times);
double st_dev(const std::vector<double>& times, double mean_time);
void test_power(int num_runs = 1000);

}