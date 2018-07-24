# Create a library called "tutrhotask" 

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${BASE_INCLUDE_DIRECTORIES} 
${CMAKE_SOURCE_DIR}/field 
#${CMAKE_SOURCE_DIR}/trackbase 
${CMAKE_SOURCE_DIR}/pnddata
${CMAKE_SOURCE_DIR}/pnddata/PidData
${CMAKE_SOURCE_DIR}/passive
${CMAKE_SOURCE_DIR}/tpc
${CMAKE_SOURCE_DIR}/lhetrack
${CMAKE_SOURCE_DIR}/emc/EmcData
${CMAKE_SOURCE_DIR}/emc/EmcReco
${CMAKE_SOURCE_DIR}/fsim
${CMAKE_SOURCE_DIR}/rho
${CMAKE_SOURCE_DIR}/rho/RhoBase
${CMAKE_SOURCE_DIR}/rho/RhoSelector
${CMAKE_SOURCE_DIR}/rho/RhoFitter
${CMAKE_SOURCE_DIR}/tutorials/rho
${CMAKE_SOURCE_DIR}/PndTools/AnalysisTools
${CMAKE_SOURCE_DIR}/PndTools/PndParticleFitters
${CMAKE_SOURCE_DIR}/genfit
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS
PndTutAnaTask.cxx
)

set(LINKDEF tutrhoLinkDef.h)
set(LIBRARY_NAME tutrhotask)

set(DEPENDENCIES PndData Rho AnalysisTools)

GENERATE_LIBRARY()
