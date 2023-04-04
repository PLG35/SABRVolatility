#-------------------------------------------------
#
# Project created by QtCreator 2018-11-01T20:51:55
#
#-------------------------------------------------

QT       += core printsupport

TARGET = launch_tests_SABR
TEMPLATE = app

SOURCES += \
    test_smile_curve.cpp \
    lib_smile_curve.cpp \
    main.cpp \
    test_planner_v2.cpp \
    test_planner_v1.cpp \
    lib_calibration.cpp \
    test_planner_calibrator.cpp

HEADERS  += test_smile_curve.h \
    lib_smile_curve.h \
    test_planner_v2.h \
    test_planner_v1.h \
    lib_calibration.h \
    test_planner_calibrator.h
