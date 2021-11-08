#include "SZLOrderedPartitionedZoneWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <new>
#include <sstream>
#include <utility>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/ref.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/unordered_set.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "FieldData.h"
#include "ItemAddress.h"
#include "ItemSetIterator.h"
#include "SZLOrderedPartitionWriter.h"
#include "writeValueArray.h"
#include "ZoneInfoCache.h"
namespace tecplot { namespace ___3933 { SZLOrderedPartitionedZoneWriter::SZLOrderedPartitionedZoneWriter( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36, ZoneInfoCache&                zoneInfoCache) : ___4709(varIter, zone, ___341, ___4564, ___4499, ___36) , m_headerWriter( varIter, zone, ___341, ___36, m_partitionFileNums, m_partitionHeaderFilePositions, m_partitionMinNodeNumbers, m_partitionMaxNodeNumbers, m_varPartitionMinMaxes) , ___2680(zoneInfoCache) , m_partitionTecUtil(___36, zone + 1) { REQUIRE(0 <= zone && ___36.___4638(zone + 1)); REQUIRE(VALID_BOOLEAN(___4499)); REQUIRE(___36.zoneIsPartitioned(zone + 1)); size_t const numVarsToWrite = static_cast<size_t>(m_varIter.___2812()); size_t const numPartitions = static_cast<size_t>(___36.zoneGetNumPartitions(zone + 1)); if (!m_partitionFileNums.alloc(numPartitions, 0) || !m_partitionHeaderFilePositions.alloc(numPartitions, ___330) || !m_partitionMinNodeNumbers.alloc(numPartitions, 0) || !m_partitionMaxNodeNumbers.alloc(numPartitions, 0) || !___3356(m_varPartitionMinMaxes, numVarsToWrite, numPartitions)) throw std::bad_alloc(); } SZLOrderedPartitionedZoneWriter::~SZLOrderedPartitionedZoneWriter() {} void SZLOrderedPartitionedZoneWriter::getCellMinMaxes( std::vector<___2479>& cellMinMaxes, ___2227 ___462, ___1844 const& dimensions, std::vector<___1352> const& fieldDatas) { REQUIRE(cellMinMaxes.size() == fieldDatas.size()); for (size_t ___4336 = 0; ___4336 < fieldDatas.size(); ++___4336) cellMinMaxes[___4336].invalidate(); for (int cellI = 0; cellI <= 1; ++cellI) { for (int cellJ = 0; cellJ <= 1; ++cellJ) { for (int cellK = 0; cellK <= 1; ++cellK) { ___2227 ___2716 = ___462 + (cellK * dimensions.___2105() + cellJ) * dimensions.i() + cellI + 1; for (size_t ___4336 = 0; ___4336 < fieldDatas.size(); ++___4336) { double ___4298 = 0.0; if (fieldDatas[___4336].___2067()) ___4298 = fieldDatas[___4336].___1780(___2716); cellMinMaxes[___4336].include(___4298); } } } } } namespace{ void applyCellMinMaxToNeighborNodeSubzones( ___2227 i, ___2227 ___2105, ___2227 ___2134, std::vector<___2479> const& cellMinMaxes, ___1844 const& partitionOffsetIJK, ___1844 const& neighborOffsetIJK, ___1881& neighborInfo) { for (___2227 ___2158 = 0; ___2158 <= 1; ++___2158) { ___2227 neighborK = ___2158 + ___2134 + partitionOffsetIJK.___2134() - neighborOffsetIJK.___2134(); if (0 <= neighborK && neighborK < neighborInfo.___2714().___2134()) { for (___2227 ___2113 = 0; ___2113 <= 1; ++___2113) { ___2227 neighborJ = ___2113 + ___2105 + partitionOffsetIJK.___2105() - neighborOffsetIJK.___2105(); if (0 <= neighborJ && neighborJ < neighborInfo.___2714().___2105()) { for (___2227 ___1841 = 0; ___1841 <= 1; ++___1841) { ___2227 neighborI = ___1841 + i + partitionOffsetIJK.i() - neighborOffsetIJK.i(); if (0 <= neighborI && neighborI < neighborInfo.___2714().i()) { ___1844 neighborNodeIJK((___81)neighborI, (___81)neighborJ, (___81)neighborK); ___2090::SubzoneOffset_t ___2734 = neighborInfo.nszAtNodeIJK(neighborNodeIJK).subzoneOffset(); neighborInfo.includeNszVarMinMax(___2734, cellMinMaxes); } } } } } } } } void SZLOrderedPartitionedZoneWriter::applyCellMinMaxesToNeighborsInRange( std::vector<___1864> const& neighborItems, ___2090::___2980 ___2977, ___1855 const& partitionRange, ___1844 const& partitionOffsetIJK, ___1844 const& partitionDimensionsIJK, std::vector<___1352> const& fieldDatas, std::vector<boost::shared_ptr<___1881> >& partitionInfos) { BOOST_FOREACH(___1864 const& neighborItem, neighborItems) { ___2090::___2980 neighborPartition = neighborItem.second; if (neighborPartition == ___2977) continue; ___1855 intersectionRange; ___1855 neighborRange = neighborItem.first; boost::geometry::intersection(partitionRange, neighborRange, intersectionRange); throwIfBadIntersectionRange(intersectionRange, ___2977, neighborPartition); ___2227 cellIMin = std::max((___2227)0, (___2227)(intersectionRange.min_corner().get<0>() - partitionOffsetIJK.i() - 1)); ___2227 ___461 = std::min((___2227)(partitionDimensionsIJK.i() - 2), (___2227)(intersectionRange.max_corner().get<0>() - partitionOffsetIJK.i())); ___2227 cellJMin = std::max((___2227)0, (___2227)(intersectionRange.min_corner().get<1>() - partitionOffsetIJK.___2105() - 1)); ___2227 ___466 = std::min((___2227)(partitionDimensionsIJK.___2105() - 2), (___2227)(intersectionRange.max_corner().get<1>() - partitionOffsetIJK.___2105())); ___2227 cellKMin = std::max((___2227)0, intersectionRange.min_corner().get<2>() - partitionOffsetIJK.___2134() - 1); ___2227 ___467 = std::min((___2227)(partitionDimensionsIJK.___2134() - 2), (___2227)(intersectionRange.max_corner().get<2>() - partitionOffsetIJK.___2134())); ___1844 neighborOffsetIJK((___81)neighborRange.min_corner().get<0>(), (___81)neighborRange.min_corner().get<1>(), (___81)neighborRange.min_corner().get<2>());
for (___2227 ___2134 = cellKMin; ___2134 <= ___467; ++___2134) { for (___2227 ___2105 = cellJMin; ___2105 <= ___466; ++___2105) { for (___2227 i = cellIMin; i <= ___461; ++i) { std::vector<___2479> cellMinMaxes(fieldDatas.size()); ___2227 ___462 = (___2134 * partitionDimensionsIJK.___2105() + ___2105) * partitionDimensionsIJK.i() + i; getCellMinMaxes(cellMinMaxes, ___462, partitionDimensionsIJK, fieldDatas); applyCellMinMaxToNeighborNodeSubzones(i, ___2105, ___2134, cellMinMaxes, partitionOffsetIJK, neighborOffsetIJK, *partitionInfos[neighborPartition]); } } } } } void SZLOrderedPartitionedZoneWriter::applyCellMinMaxesToNeighborNodeSubzones( ___2090::___2980 ___2977, std::vector<___1352> const& nodalFieldDatas, std::vector<boost::shared_ptr<___1881> >& partitionInfos, ___1863 const& partitionTree) { REQUIRE(___2977 >= 0); ___1844 partitionOffsetIJK; ___2337.zonePartitionGetIJKOffset(___2677 + 1, ___2977 + 1, partitionOffsetIJK); ___1853 ___2474(partitionOffsetIJK.i(), partitionOffsetIJK.___2105(), partitionOffsetIJK.___2134()); ___1844 partitionDimensionsIJK; ___2337.zonePartitionGetIJK(___2677 + 1, ___2977 + 1, partitionDimensionsIJK); ___1844 partitionMaxIJK = partitionOffsetIJK + partitionDimensionsIJK - 1; ___1853 ___2364(partitionMaxIJK.i(), partitionMaxIJK.___2105(), partitionMaxIJK.___2134()); std::vector<___1864> neighborItems; ___1853 iFaceMinCorner(partitionMaxIJK.i(), partitionOffsetIJK.___2105(), partitionOffsetIJK.___2134()); ___1855 iFaceRange(iFaceMinCorner, ___2364); partitionTree.query(boost::geometry::index::intersects(iFaceRange), std::back_inserter(neighborItems)); applyCellMinMaxesToNeighborsInRange(neighborItems, ___2977, iFaceRange, partitionOffsetIJK, partitionDimensionsIJK, nodalFieldDatas, partitionInfos); neighborItems.clear(); ___1853 jFaceMinCorner(partitionOffsetIJK.i(), partitionMaxIJK.___2105(), partitionOffsetIJK.___2134()); ___1855 jFaceRange(jFaceMinCorner, ___2364); partitionTree.query(boost::geometry::index::intersects(jFaceRange), std::back_inserter(neighborItems)); applyCellMinMaxesToNeighborsInRange(neighborItems, ___2977, jFaceRange, partitionOffsetIJK, partitionDimensionsIJK, nodalFieldDatas, partitionInfos); neighborItems.clear(); ___1853 kFaceMinCorner(partitionOffsetIJK.i(), partitionOffsetIJK.___2105(), partitionMaxIJK.___2134()); ___1855 kFaceRange(kFaceMinCorner, ___2364); partitionTree.query(boost::geometry::index::intersects(kFaceRange), std::back_inserter(neighborItems)); applyCellMinMaxesToNeighborsInRange(neighborItems, ___2977, kFaceRange, partitionOffsetIJK, partitionDimensionsIJK, nodalFieldDatas, partitionInfos); } void SZLOrderedPartitionedZoneWriter::retrieveNodalFieldDataPtrsForPartition( ___37& partitionTecUtilDecorator, ___2090::___2980 ___2977, std::vector<___1352> &nodalFieldDatas) { m_varIter.reset(); ___4352 const baseVar = m_varIter.baseItem(); while (m_varIter.hasNext()) { ___4352 const datasetVar = m_varIter.next(); ___4352 const fileVar = datasetVar - baseVar; if (___2337.___4353(datasetVar + 1) && !___2337.___926(___2677 + 1, datasetVar + 1)) { if (___2337.___910(___2677 + 1, datasetVar + 1) == ___4330) nodalFieldDatas[fileVar] = ___1352(&partitionTecUtilDecorator, ___2977 + 1, datasetVar + 1, false, false); else nodalFieldDatas[fileVar] = ___1352(&partitionTecUtilDecorator, ___2977 + 1, datasetVar + 1, false, true); } } } void SZLOrderedPartitionedZoneWriter::throwIfBadIntersectionRange( ___1855 const& intersectionRange, ___2090::___2980 ___2977, ___2090::___2980 neighborPartition) { if (intersectionRange.min_corner().get<0>() != intersectionRange.max_corner().get<0>() && intersectionRange.min_corner().get<1>() != intersectionRange.max_corner().get<1>() && intersectionRange.min_corner().get<2>() != intersectionRange.max_corner().get<2>()) { std::ostringstream ___2892; ___2892 << "Error writing zone " << ___2677 + 1 << ": partition " << ___2977 + 1 << " overlaps partition " << neighborPartition << " by more than one node layer. The overlap is (" << intersectionRange.min_corner().get<0>() << '-' << intersectionRange.max_corner().get<0>() << ", " << intersectionRange.min_corner().get<1>() << '-' << intersectionRange.max_corner().get<1>() << ", " << intersectionRange.min_corner().get<2>() << '-' << intersectionRange.max_corner().get<2>() << ")." << " Ghost cells cannot be output for ordered zones. Please correct the partition index ranges in your calls to tecijkptn."; throw std::runtime_error(___2892.str()); } } void SZLOrderedPartitionedZoneWriter::exchangeGhostInfo( std::vector<boost::shared_ptr<___1881> >& partitionInfos, std::vector<___1864> const& ___2981) { ___1863 partitionTree(___2981); PartitionTecUtilDecorator partitionTecUtilDecorator(___2337, ___2677 + 1); ___4352 const numVarsToWrite = m_varIter.___2812(); ___2090::___2980 const numPartitions = static_cast<___2090::___2980>(___2337.zoneGetNumPartitions(___2677 + 1));
std::vector<___1864> neighborPartitions; for (___2090::___2980 ___2977 = 0; ___2977 < numPartitions; ++___2977) { std::vector<___1352> nodalFieldDatas(numVarsToWrite); retrieveNodalFieldDataPtrsForPartition(partitionTecUtilDecorator, ___2977, nodalFieldDatas); applyCellMinMaxesToNeighborNodeSubzones(___2977, nodalFieldDatas, partitionInfos, partitionTree); } } void SZLOrderedPartitionedZoneWriter::getPartitionExtentsWithGhostNodes( ___2090::___2980 ___2977, ___1844& partitionMinIJK, ___1844& partitionMaxIJK) { REQUIRE(0 <= ___2977 && ___2977 < static_cast<___2090::___2980>(___2337.zoneGetNumPartitions(___2677 + 1))); ___1844 partitionSize; ___2337.zonePartitionGetIJK(___2677 + 1, ___2977 + 1, partitionSize); ___2337.zonePartitionGetIJKOffset(___2677 + 1, ___2977 + 1, partitionMinIJK); partitionMaxIJK = partitionMinIJK + partitionSize - 1; } void SZLOrderedPartitionedZoneWriter::trimGhostNodes(___1844 &partitionMaxIJK) { ___1844 zoneSize; ___2337.___4615(___2677 + 1, zoneSize); if (partitionMaxIJK.i() < zoneSize.i() - 1) partitionMaxIJK.setI(partitionMaxIJK.i() - 1); if (partitionMaxIJK.___2105() < zoneSize.___2105() - 1) partitionMaxIJK.setJ(partitionMaxIJK.___2105() - 1); if (partitionMaxIJK.___2134() < zoneSize.___2134() - 1) partitionMaxIJK.___3497(partitionMaxIJK.___2134() - 1); } void SZLOrderedPartitionedZoneWriter::getPartitionExtentsWithoutGhostNodes( ___2090::___2980 ___2977, ___1844& partitionMinIJK, ___1844& partitionMaxIJK) { REQUIRE(0 <= ___2977 && ___2977 < static_cast<___2090::___2980>(___2337.zoneGetNumPartitions(___2677 + 1))); getPartitionExtentsWithGhostNodes(___2977, partitionMinIJK, partitionMaxIJK); trimGhostNodes(partitionMaxIJK); } void SZLOrderedPartitionedZoneWriter::createPartitionWriters() { std::vector<boost::shared_ptr<___1881> > partitionInfos; ___2090::___2980 numPartitions = static_cast<___2090::___2980>(___2337.zoneGetNumPartitions(___2677 + 1)); std::vector<___1864> ___2981; ___2981.reserve((size_t)numPartitions); ___1844 zoneSize; ___2337.___4615(___2677 + 1, zoneSize); for (___2090::___2980 ___2977 = 0; ___2977 < numPartitions; ++___2977) { partitionInfos.push_back(___2680.getIJKZonePartitionInfo(___2677, m_baseZone, ___2977)); ___1844 partitionMinIJK; ___1844 partitionMaxIJK; getPartitionExtentsWithoutGhostNodes(___2977, partitionMinIJK, partitionMaxIJK); m_partitionMinNodeNumbers[___2977] = zoneSize.offsetAtIJK(partitionMinIJK); m_partitionMaxNodeNumbers[___2977] = zoneSize.offsetAtIJK(partitionMaxIJK); ___1853 ___2478(partitionMinIJK.i(), partitionMinIJK.___2105(), partitionMinIJK.___2134()); ___1853 ___2372(partitionMaxIJK.i(), partitionMaxIJK.___2105(), partitionMaxIJK.___2134()); ___2981.push_back(std::make_pair(___1855(___2478, ___2372), ___2977)); } exchangeGhostInfo(partitionInfos, ___2981); for (___2090::___2980 ___2977 = 0; ___2977 < numPartitions; ++___2977) { m_partitionWriters[___2977] = boost::make_shared<SZLOrderedPartitionWriter> ( boost::ref(m_varIter), ___2677, m_baseZone, ___2977, boost::ref(m_writeVariables), m_writeConnectivity, boost::ref(m_partitionTecUtil), partitionInfos[___2977]); } } ___372 SZLOrderedPartitionedZoneWriter::writeZoneData(FileWriterInterface& szpltFile) { if (m_partitionWriters.empty()) createPartitionWriters(); for(___4636 ___2977 = 0; ___2977 < ___2337.zoneGetNumPartitions(___2677 + 1); ++___2977) { m_partitionWriters[___2977]->writeZone(szpltFile, szpltFile.fileLoc()); m_partitionHeaderFilePositions[___2977] = m_partitionWriters[___2977]->getZoneHeaderFilePosition(); m_varIter.reset(); ___4352 const baseVar = m_varIter.baseItem(); while (m_varIter.hasNext()) { ___4352 const datasetVar = m_varIter.next(); ___4352 const fileVar = datasetVar - baseVar; m_varPartitionMinMaxes[fileVar][___2977] = m_partitionWriters[___2977]->varMinMax(datasetVar); } } return ___4226; } uint64_t SZLOrderedPartitionedZoneWriter::zoneDataFileSize(bool ___2002) { if (m_partitionWriters.empty()) createPartitionWriters(); uint64_t ___3358 = 0; for(___4636 ___2977 = 0; ___2977 < ___2337.zoneGetNumPartitions(___2677 + 1); ++___2977) ___3358 += m_partitionWriters[___2977]->zoneFileSize(___2002); return ___3358; } ___372 SZLOrderedPartitionedZoneWriter::writeZoneConnectivity(FileWriterInterface&  ) { return ___4226; } uint64_t SZLOrderedPartitionedZoneWriter::zoneConnectivityFileSize(bool  ) { return 0; } ___372 SZLOrderedPartitionedZoneWriter::writeZoneHeader(FileWriterInterface& szpltFile) { REQUIRE(szpltFile.___2041()); return m_headerWriter.write(szpltFile); } uint64_t SZLOrderedPartitionedZoneWriter::zoneHeaderFileSize(bool ___2002) { return m_headerWriter.sizeInFile(___2002); } }}