/*
 * tests.h
 *
 *  Created on: 2 ���. 2021 �.
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
#include <thread>
#include <chrono>
#include "gwell.h"
#include <fstream>

std::ostream& operator<<(std::ostream& os, const PointXYZV& p);

void TestCinco();

void TestPdLaplParallel();

void ShowPdXY();

