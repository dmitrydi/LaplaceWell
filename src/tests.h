/*
 * tests.h
 *
 *  Created on: 2 мар. 2021 г.
 *      Author: Dmitry_Di
 */

#pragma once

#include "auxillary.h"
#include "test_runner.h"
#include "profile.h"
#include <vector>
#include <iostream>
#include <limits>
#include <iterator>
#include <cmath>
#include <future>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>
#include "gwell.h"
#include <fstream>

std::ostream& operator<<(std::ostream& os, const PointXYZV& p);

void TestMakeGrid();

void TestNPaginate();
void TestNPaginateXY(size_t nx, size_t ny, size_t nthread);

MatrixXYZV _TestParallelDum(const std::vector<double>& xs, const std::vector<double>& ys, int nthread);
void TestParallelDum();

void TestPDXY();

void TestMatrix3DV();

void TestMatrix3DV_Pages();

void Test_pd_lapl_m();

void TestPDXY_old_new();

void TestCincoNew();

