target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/RankZeroTensorOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RankOneTensorOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RankTwoTensorOperator.cpp
    )

add_subdirectory(one_electron)
add_subdirectory(two_electron)
