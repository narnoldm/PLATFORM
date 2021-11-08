#include "SzlFileLoader.h"
#include "ThirdPartyHeadersBegin.h"
#include <algorithm>
#include <math.h>
#include "ThirdPartyHeadersEnd.h"
#include "FileDescription.h"
#include "ORBFESubzonePartitioner.h"
#include "AltTecUtil.h"
#include "stringformat.h"
#include "zoneUtil.h"
namespace tecplot { namespace ___3933 { OrbFESubzonePartitioner::OrbFESubzonePartitioner( ___37&               ___36, ___4636               zone, ___2090::ItemOffset_t fixedSubzoneSize, bool                      sortItems) : ___2677(zone) , m_fixedSubzoneSize(fixedSubzoneSize) , m_cellOrb(zone, OrthogonalBisection::BisectionType_ZoneCells, fixedSubzoneSize, sortItems) , m_nodeOrb(zone, OrthogonalBisection::BisectionType_ZoneNodes, fixedSubzoneSize, sortItems) { REQUIRE(___36.___4638(___2677 + 1)); REQUIRE(m_fixedSubzoneSize <= ___2090::MAX_ITEM_OFFSET+1); REQUIRE(___3894(___36.___4620(___2677 + 1))); partitionIntoSubzones(___36); } OrbFESubzonePartitioner::~OrbFESubzonePartitioner() { m_nodeOrb.___937(); m_cellOrb.___937(); m_szCoordsOfOrginalZoneCells.___937(); m_szCoordsOfOrginalZoneNodes.___937(); } ___372 OrbFESubzonePartitioner::partitionIntoSubzones(___37& ___36) { bool ___2039 = true; size_t const messageSize = 200; char statusMessage[messageSize]; snprintf(statusMessage, messageSize, "Determining node subzones for zone %" PRIu64 "...", uint64_t(___2677)+1); ___36.___3778(statusMessage);
 #ifdef TIME_FE_DECOMPOSITION
uint64_t startTimeInMS = ___717();
 #endif
___2039 = ___2039 && m_nodeOrb.performBisection(___36);
 #ifdef TIME_FE_DECOMPOSITION
uint64_t endTimeInMS = ___717(); ___1931("Time to create node subzones is %" PRIu64 " ms", long(endTimeInMS-startTimeInMS));
 #endif
___2039 = ___2039 && m_nodeOrb.getSzCoordByOriginalItemArray(___36, m_szCoordsOfOrginalZoneNodes); if ( ___2039 ) { size_t const messageSize = 200; char statusMessage[messageSize]; snprintf(statusMessage, messageSize, "Determining cell subzones for zone %" PRIu64 "...", uint64_t(___2677)+1); ___36.___3778(statusMessage);
 #ifdef TIME_FE_DECOMPOSITION
uint64_t startTimeInMS = ___717();
 #endif
___2039 = ___2039 && m_cellOrb.performBisection(___36);
 #ifdef TIME_FE_DECOMPOSITION
uint64_t endTimeInMS = ___717(); ___1931("Time to create cell subzones is %" PRIu64 " ms", endTimeInMS-startTimeInMS);
 #endif
___2039 = ___2039 && m_cellOrb.getSzCoordByOriginalItemArray(___36, m_szCoordsOfOrginalZoneCells); } return ___372(___2039); } ___465 OrbFESubzonePartitioner::numCellsInZone() const { ENSURE(m_cellOrb.queryNumItems() > 0); return static_cast<___465>(m_cellOrb.queryNumItems() + m_cellOrb.queryNumGhostItems()); } ___2090::SubzoneOffset_t OrbFESubzonePartitioner::___2783() const { ENSURE(m_cellOrb.queryNumberDomains() > 0); return m_cellOrb.queryNumberDomains(); } ___2090::ItemOffset_t OrbFESubzonePartitioner::___2782(___2090::SubzoneOffset_t ___469) const { ___2090::ItemOffset_t const cszSize = m_cellOrb.getDomainSize(___469); ENSURE(cszSize > 0 && cszSize <= ___2090::MAX_ITEM_OFFSET+1); return cszSize; } ___2090 OrbFESubzonePartitioner::szCoordinateAtZoneCell(___465 zoneCell) const { REQUIRE(0 <= zoneCell); ___2090 const szCoordinate = m_szCoordsOfOrginalZoneCells[zoneCell]; ENSURE(IMPLICATION(static_cast<___4636>(szCoordinate.___2977()) == ___2677, szCoordinate.subzoneOffset() < ___2783())); ENSURE(IMPLICATION(static_cast<___4636>(szCoordinate.___2977()) == ___2677, szCoordinate.itemOffset() < ___2782(szCoordinate.subzoneOffset()))); ENSURE(IMPLICATION(static_cast<___4636>(szCoordinate.___2977()) == ___2677, ___4608(szCoordinate) == zoneCell)); return szCoordinate; } ___465 OrbFESubzonePartitioner::___4608(___2090 ___451) const { REQUIRE(___451.subzoneOffset()<___2783()); REQUIRE(___451.itemOffset()<___2782(___451.subzoneOffset())); ___465 newZoneCell = ___451.subzoneOffset() * m_fixedSubzoneSize + ___451.itemOffset(); ___465 zoneCell = m_cellOrb.queryPositionbyOffset(newZoneCell); ENSURE(zoneCell < numCellsInZone()); return zoneCell; } ___2718 OrbFESubzonePartitioner::numNodesInZone() const { ENSURE(m_nodeOrb.queryNumItems() > 0); return static_cast<___2718>(m_nodeOrb.queryNumItems() + m_nodeOrb.queryNumGhostItems()); } ___2090::SubzoneOffset_t OrbFESubzonePartitioner::___2823() const { ENSURE(m_nodeOrb.queryNumberDomains() > 0); return m_nodeOrb.queryNumberDomains(); } ___2090::ItemOffset_t OrbFESubzonePartitioner::___2822(___2090::SubzoneOffset_t ___2734) const { ___2090::ItemOffset_t const nszSize = m_nodeOrb.getDomainSize(___2734); ENSURE(nszSize > 0 && nszSize <= ___2090::MAX_ITEM_OFFSET+1); return nszSize; } ___2090 OrbFESubzonePartitioner::___3924(___2718 ___4656) const { ___2090 const nodeAddress = m_szCoordsOfOrginalZoneNodes[___4656]; ENSURE(IMPLICATION(static_cast<___4636>(nodeAddress.___2977()) == ___2677, nodeAddress.subzoneOffset() < ___2823())); ENSURE(IMPLICATION(static_cast<___4636>(nodeAddress.___2977()) == ___2677, nodeAddress.itemOffset() < ___2822(nodeAddress.subzoneOffset()))); ENSURE(IMPLICATION(static_cast<___4636>(nodeAddress.___2977()) == ___2677, ___4657(nodeAddress) == ___4656)); return nodeAddress; } ___2718 OrbFESubzonePartitioner::___4657(___2090 nodeAddress) const { REQUIRE(nodeAddress.subzoneOffset()<___2823()); REQUIRE(nodeAddress.itemOffset()<___2822(nodeAddress.subzoneOffset())); ___2718 const newZoneNode = nodeAddress.subzoneOffset() * m_fixedSubzoneSize + nodeAddress.itemOffset(); ___2718 const ___4656 =  ___2718( m_nodeOrb.queryPositionbyOffset(newZoneNode) ); ENSURE(___4656 < numNodesInZone()); return ___4656; } void OrbFESubzonePartitioner::setNodeSubzoneCoordinate(___2718 ___4656, ___2090 ___2759) { REQUIRE(0 <= ___4656 && ___4656 < numNodesInZone()); REQUIRE(___2759.___14() == ___2090::SzlAddressType); m_szCoordsOfOrginalZoneNodes[___4656] = ___2759; } }}