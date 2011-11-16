CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(tmv)

ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(test EXCLUDE_FROM_ALL)

ADD_SUBDIRECTORY(examples EXCLUDE_FROM_ALL)

INSTALL(DIRECTORY include/ DESTINATION include)

ADD_CUSTOM_TARGET(test 
    DEPENDS tmvtest1 tmvtest2 tmvtest3
    WORKING_DIRECTORY test)

ADD_CUSTOM_TARGET(smalltest 
    DEPENDS tmvtest1a tmvtest1b tmvtest1c tmvtest1d tmvtest2a tmvtest2b tmvtest2c tmvtest3a tmvtest3b tmvtest3c tmvtest3d tmvtest3e
    WORKING_DIRECTORY test)

ADD_CUSTOM_TARGET(examples
    DEPENDS vector matrix division bandmatrix symmatrix
    WORKING_DIRECTORY examples)

