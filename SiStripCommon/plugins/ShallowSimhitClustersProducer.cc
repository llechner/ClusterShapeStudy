#include "CalibTracker/SiStripCommon/interface/ShallowSimhitClustersProducer.h"

#include "CalibTracker/SiStripCommon/interface/ShallowTools.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "CondFormats/SiStripObjects/interface/SiStripLorentzAngle.h"
#include "CondFormats/DataRecord/interface/SiStripLorentzAngleRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

ShallowSimhitClustersProducer::ShallowSimhitClustersProducer(const edm::ParameterSet& iConfig)
  : clusters_token_( consumes< edmNew::DetSetVector<SiStripCluster> >(iConfig.getParameter<edm::InputTag>("Clusters"))),
    Prefix( iConfig.getParameter<std::string>("Prefix") ),
		runningmode_( iConfig.getParameter<std::string>("runningMode") )
{
	std::vector<edm::InputTag> simhits_tags = iConfig.getParameter<std::vector<edm::InputTag> >("InputTags");
	for(auto itag : simhits_tags) {
		simhits_tokens_.push_back(
			consumes< std::vector<PSimHit> >(itag)
			);
	}

  produces <std::vector<unsigned int> >  ( Prefix + "trackMother"      );
  produces <std::vector<unsigned> >      ( Prefix + "hits"       );
  produces <std::vector<float> >         ( Prefix + "strip"      );
  produces <std::vector<float> >         ( Prefix + "localtheta" );
  produces <std::vector<float> >         ( Prefix + "localphi"   );
  produces <std::vector<float> >         ( Prefix + "localx"     );
  produces <std::vector<float> >         ( Prefix + "localy"     );
  produces <std::vector<float> >         ( Prefix + "localz"     );
  produces <std::vector<float> >         ( Prefix + "globalx"     );
  produces <std::vector<float> >         ( Prefix + "globaly"     );
  produces <std::vector<float> >         ( Prefix + "globalz"     );
  produces <std::vector<float> >         ( Prefix + "momentum"   );
  produces <std::vector<float> >         ( Prefix + "energyloss" );
  produces <std::vector<float> >         ( Prefix + "time"       );
  produces <std::vector<int> >            ( Prefix + "particle"   );
  produces <std::vector<unsigned short> > ( Prefix + "process"    );

  produces <std::vector<unsigned int> >  ( Prefix + "trackMother2"      );
  produces <std::vector<float> >         ( Prefix + "strip2"      );
  produces <std::vector<float> >         ( Prefix + "localtheta2" );
  produces <std::vector<float> >         ( Prefix + "localphi2"   );
  produces <std::vector<float> >         ( Prefix + "localx2"     );
  produces <std::vector<float> >         ( Prefix + "localy2"     );
  produces <std::vector<float> >         ( Prefix + "localz2"     );
  produces <std::vector<float> >         ( Prefix + "globalx2"     );
  produces <std::vector<float> >         ( Prefix + "globaly2"     );
  produces <std::vector<float> >         ( Prefix + "globalz2"     );
  produces <std::vector<float> >         ( Prefix + "momentum2"   );
  produces <std::vector<float> >         ( Prefix + "energyloss2" );
  produces <std::vector<float> >         ( Prefix + "time2"       );
  produces <std::vector<int> >            ( Prefix + "particle2"   );
  produces <std::vector<unsigned short> > ( Prefix + "process2"    );

  produces <std::vector<unsigned int> >  ( Prefix + "trackMother3"      );
  produces <std::vector<float> >         ( Prefix + "localx3"     );
  produces <std::vector<float> >         ( Prefix + "localy3"     );
  produces <std::vector<float> >         ( Prefix + "localz3"     );
  produces <std::vector<float> >         ( Prefix + "energyloss3" );
  produces <std::vector<float> >         ( Prefix + "time3"       );
  produces <std::vector<int> >           ( Prefix + "particle3"   );

  produces <std::vector<unsigned int> >  ( Prefix + "trackMother4"      );
  produces <std::vector<float> >         ( Prefix + "localx4"     );
  produces <std::vector<float> >         ( Prefix + "localy4"     );
  produces <std::vector<float> >         ( Prefix + "localz4"     );
  produces <std::vector<float> >         ( Prefix + "energyloss4" );
  produces <std::vector<float> >         ( Prefix + "time4"       );
  produces <std::vector<int> >           ( Prefix + "particle4"   );

}	      

void ShallowSimhitClustersProducer::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  shallow::CLUSTERMAP clustermap = shallow::make_cluster_map(iEvent,clusters_token_);

  int size = clustermap.size();
  auto          trackMother  = std::make_unique <std::vector<unsigned int>>       (size, 0);
  auto          hits         = std::make_unique <std::vector<unsigned>>       (size, 0);
  auto          strip        = std::make_unique <std::vector<float>>          (size, -100);
  auto          localtheta   = std::make_unique <std::vector<float>>          (size, -100);
  auto          localphi     = std::make_unique <std::vector<float>>          (size, -100);
  auto          localx       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localy       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localz       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globalx       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globaly       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globalz       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          momentum     = std::make_unique <std::vector<float>>          (size, -1);
  auto          energyloss   = std::make_unique <std::vector<float>>          (size, -100);
  auto          time         = std::make_unique <std::vector<float>>          (size, -1);
  auto          particle     = std::make_unique <std::vector<int>>            (size, -500);
  auto          process      = std::make_unique <std::vector<unsigned short>> (size, 0);

  auto          trackMother2  = std::make_unique <std::vector<unsigned int>>       (size, 0);
  auto          strip2        = std::make_unique <std::vector<float>>          (size, -100);
  auto          localtheta2   = std::make_unique <std::vector<float>>          (size, -100);
  auto          localphi2     = std::make_unique <std::vector<float>>          (size, -100);
  auto          localx2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localy2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localz2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globalx2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globaly2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          globalz2       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          momentum2     = std::make_unique <std::vector<float>>          (size, -1);
  auto          energyloss2   = std::make_unique <std::vector<float>>          (size, -100);
  auto          time2         = std::make_unique <std::vector<float>>          (size, -1);
  auto          particle2     = std::make_unique <std::vector<int>>            (size, -500);
  auto          process2      = std::make_unique <std::vector<unsigned short>> (size, 0);

  auto          trackMother3  = std::make_unique <std::vector<unsigned int>>       (size, 0);
  auto          localx3       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localy3       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localz3       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          energyloss3   = std::make_unique <std::vector<float>>          (size, -100);
  auto          time3         = std::make_unique <std::vector<float>>          (size, -1);
  auto          particle3     = std::make_unique <std::vector<int>>            (size, -500);

  auto          trackMother4  = std::make_unique <std::vector<unsigned int>>       (size, 0);
  auto          localx4       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localy4       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          localz4       = std::make_unique <std::vector<float>>          (size, -1000);
  auto          energyloss4   = std::make_unique <std::vector<float>>          (size, -100);
  auto          time4         = std::make_unique <std::vector<float>>          (size, -1);
  auto          particle4     = std::make_unique <std::vector<int>>            (size, -500);

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  edm::ESHandle<MagneticField> magfield;		     iSetup.get<IdealMagneticFieldRecord>().get(magfield);		      
  edm::ESHandle<SiStripLorentzAngle> SiStripLorentzAngle;    iSetup.get<SiStripLorentzAngleRcd>().get(runningmode_, SiStripLorentzAngle);      
  edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;  iEvent.getByLabel("siStripClusters", "", clusters);

	for(auto& simhit_token : simhits_tokens_) {
		edm::Handle< std::vector<PSimHit> > simhits; 
		iEvent.getByToken(simhit_token, simhits);
    for(auto const& hit : *simhits ) {
      
      const uint32_t id = hit.detUnitId();
      const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( id ) );
      const LocalVector drift = shallow::drift(theStripDet, *magfield, *SiStripLorentzAngle);

      const float driftedstrip_ = theStripDet->specificTopology().strip( hit.localPosition()+0.5*drift );
      const float     hitstrip_ = theStripDet->specificTopology().strip( hit.localPosition()           );

      shallow::CLUSTERMAP::const_iterator cluster = match_cluster( id, driftedstrip_, clustermap, *clusters); 
      if(cluster != clustermap.end()) {
				unsigned i = cluster->second;
                if (hit.particleType() == 22) continue; // remove photons from simulated hits
				hits->at(i)+=1;
				if(hits->at(i) == 1) {
					trackMother->at(i) = hit.trackId();
					strip->at(i) = hitstrip_;
					localtheta->at(i) = hit.thetaAtEntry();
					localphi->at(i) = hit.phiAtEntry();
					localx->at(i) = hit.localPosition().x();
					localy->at(i) = hit.localPosition().y();
					localz->at(i) = hit.localPosition().z();
                    globalx->at(i)  = theStripDet->toGlobal(hit.localPosition()).x();
                    globaly->at(i)  = theStripDet->toGlobal(hit.localPosition()).y();
                    globalz->at(i)  = theStripDet->toGlobal(hit.localPosition()).z();
					momentum->at(i) = hit.pabs();
					energyloss->at(i) = hit.energyLoss();
					time->at(i) = hit.timeOfFlight();
					particle->at(i) = hit.particleType();
					process->at(i) = hit.processType();
				} else if(hits->at(i) == 2) {
					trackMother2->at(i) = hit.trackId();
					strip2->at(i) = hitstrip_;
					localtheta2->at(i) = hit.thetaAtEntry();
					localphi2->at(i) = hit.phiAtEntry();
					localx2->at(i) = hit.localPosition().x();
					localy2->at(i) = hit.localPosition().y();
					localz2->at(i) = hit.localPosition().z();
                    globalx2->at(i)  = theStripDet->toGlobal(hit.localPosition()).x();
                    globaly2->at(i)  = theStripDet->toGlobal(hit.localPosition()).y();
                    globalz2->at(i)  = theStripDet->toGlobal(hit.localPosition()).z();
					momentum2->at(i) = hit.pabs();
					energyloss2->at(i) = hit.energyLoss();
					time2->at(i) = hit.timeOfFlight();
					particle2->at(i) = hit.particleType();
					process2->at(i) = hit.processType();
				} else if(hits->at(i) == 3) {
					trackMother3->at(i) = hit.trackId();
					localx3->at(i) = hit.localPosition().x();
					localy3->at(i) = hit.localPosition().y();
					localz3->at(i) = hit.localPosition().z();
					energyloss3->at(i) = hit.energyLoss();
					time3->at(i) = hit.timeOfFlight();
					particle3->at(i) = hit.particleType();
				} else if(hits->at(i) == 4) {
					trackMother4->at(i) = hit.trackId();
					localx4->at(i) = hit.localPosition().x();
					localy4->at(i) = hit.localPosition().y();
					localz4->at(i) = hit.localPosition().z();
					energyloss4->at(i) = hit.energyLoss();
					time4->at(i) = hit.timeOfFlight();
					particle4->at(i) = hit.particleType();
				}

      }    
    } 
  }

  iEvent.put(std::move(trackMother),        Prefix + "trackMother"      );
  iEvent.put(std::move(hits),        Prefix + "hits"      );
  iEvent.put(std::move(strip),       Prefix + "strip"      );
  iEvent.put(std::move(localtheta),  Prefix + "localtheta" );
  iEvent.put(std::move(localphi),    Prefix + "localphi" );
  iEvent.put(std::move(localx),      Prefix + "localx" );
  iEvent.put(std::move(localy),      Prefix + "localy" );
  iEvent.put(std::move(localz),      Prefix + "localz" );
  iEvent.put(std::move(globalx),     Prefix + "globalx" );
  iEvent.put(std::move(globaly),     Prefix + "globaly" );
  iEvent.put(std::move(globalz),     Prefix + "globalz" );
  iEvent.put(std::move(momentum),    Prefix + "momentum" );
  iEvent.put(std::move(energyloss),  Prefix + "energyloss" );
  iEvent.put(std::move(time),        Prefix + "time" );
  iEvent.put(std::move(particle),    Prefix + "particle" );
  iEvent.put(std::move(process),     Prefix + "process" );

  iEvent.put(std::move(trackMother2),        Prefix + "trackMother2"      );
  iEvent.put(std::move(strip2),       Prefix + "strip2"      );
  iEvent.put(std::move(localtheta2),  Prefix + "localtheta2" );
  iEvent.put(std::move(localphi2),    Prefix + "localphi2" );
  iEvent.put(std::move(localx2),      Prefix + "localx2" );
  iEvent.put(std::move(localy2),      Prefix + "localy2" );
  iEvent.put(std::move(localz2),      Prefix + "localz2" );
  iEvent.put(std::move(globalx2),     Prefix + "globalx2" );
  iEvent.put(std::move(globaly2),     Prefix + "globaly2" );
  iEvent.put(std::move(globalz2),     Prefix + "globalz2" );
  iEvent.put(std::move(momentum2),    Prefix + "momentum2" );
  iEvent.put(std::move(energyloss2),  Prefix + "energyloss2" );
  iEvent.put(std::move(time2),        Prefix + "time2" );
  iEvent.put(std::move(particle2),    Prefix + "particle2" );
  iEvent.put(std::move(process2),     Prefix + "process2" );

  iEvent.put(std::move(trackMother3),        Prefix + "trackMother3"      );
  iEvent.put(std::move(localx3),      Prefix + "localx3" );
  iEvent.put(std::move(localy3),      Prefix + "localy3" );
  iEvent.put(std::move(localz3),      Prefix + "localz3" );
  iEvent.put(std::move(energyloss3),  Prefix + "energyloss3" );
  iEvent.put(std::move(time3),        Prefix + "time3" );
  iEvent.put(std::move(particle3),    Prefix + "particle3" );

  iEvent.put(std::move(trackMother4),        Prefix + "trackMother4"      );
  iEvent.put(std::move(localx4),      Prefix + "localx4" );
  iEvent.put(std::move(localy4),      Prefix + "localy4" );
  iEvent.put(std::move(localz4),      Prefix + "localz4" );
  iEvent.put(std::move(energyloss4),  Prefix + "energyloss4" );
  iEvent.put(std::move(time4),        Prefix + "time4" );
  iEvent.put(std::move(particle4),    Prefix + "particle4" );
}

shallow::CLUSTERMAP::const_iterator ShallowSimhitClustersProducer::
match_cluster( const unsigned& id, const float& strip_, const shallow::CLUSTERMAP& clustermap, const edmNew::DetSetVector<SiStripCluster>& clusters) const {
  shallow::CLUSTERMAP::const_iterator cluster = clustermap.end();
  edmNew::DetSetVector<SiStripCluster>::const_iterator clustersDetSet = clusters.find(id);
  if( clustersDetSet != clusters.end() ) {
    edmNew::DetSet<SiStripCluster>::const_iterator left, right=clustersDetSet->begin();
    while( right != clustersDetSet->end() && strip_ > right->barycenter() ) 
      right++;
    left = right-1;
    if(right!=clustersDetSet->end() && right!=clustersDetSet->begin()) { 
      unsigned firstStrip = (right->barycenter()-strip_) < (strip_-left->barycenter()) ? right->firstStrip() : left->firstStrip();
      cluster = clustermap.find( std::make_pair( id, firstStrip));
    }
    else if(right != clustersDetSet->begin())
      cluster = clustermap.find( std::make_pair( id, left->firstStrip()));
    else 
      cluster = clustermap.find( std::make_pair( id, right->firstStrip()));
  }
  return cluster;
}

