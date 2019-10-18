#include "CalibTracker/SiStripCommon/interface/ShallowGainCalibration.h"
#include "CalibFormats/SiStripObjects/interface/SiStripGain.h" 
#include "CalibTracker/Records/interface/SiStripGainRcd.h"  

using namespace edm;
using namespace reco;
using namespace std;

ShallowGainCalibration::ShallowGainCalibration(const edm::ParameterSet& iConfig)
  :  tracks_token_( consumes< edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("Tracks")) ),
	 association_token_( consumes< TrajTrackAssociationCollection >(iConfig.getParameter<edm::InputTag>("Tracks")) ),
     clusters_token_( consumes< edmNew::DetSetVector<SiStripCluster> >( iConfig.getParameter<edm::InputTag>("Clusters") ) ),
     Suffix       ( iConfig.getParameter<std::string>("Suffix")    ),
     Prefix       ( iConfig.getParameter<std::string>("Prefix") )
{
  produces <std::vector<int> >            ( Prefix + "trackindex"     + Suffix );
  produces <std::vector<unsigned int> >   ( Prefix + "rawid"          + Suffix );
  produces <std::vector<unsigned int> >   ( Prefix + "subdetid"       + Suffix );
  produces <std::vector<unsigned int> >   ( Prefix + "clusteridx"     + Suffix );
  produces <std::vector<unsigned int> >   ( Prefix + "stripidx"       + Suffix );
  produces <std::vector<double> >         ( Prefix + "localtheta"     + Suffix );
  produces <std::vector<double> >         ( Prefix + "localphi"       + Suffix );
  produces <std::vector<double> >         ( Prefix + "localdirx"      + Suffix );
  produces <std::vector<double> >         ( Prefix + "localdiry"      + Suffix );
  produces <std::vector<double> >         ( Prefix + "localdirz"      + Suffix );
  produces <std::vector<unsigned short> > ( Prefix + "firststrip"     + Suffix );
  produces <std::vector<unsigned short> > ( Prefix + "nstrips"        + Suffix );
  produces <std::vector<bool> >           ( Prefix + "saturation"     + Suffix );
  produces <std::vector<bool> >           ( Prefix + "overlapping"    + Suffix );
  produces <std::vector<bool> >           ( Prefix + "farfromedge"    + Suffix );
  produces <std::vector<int> >            ( Prefix + "layerwheel"     + Suffix );
  produces <std::vector<unsigned int> >   ( Prefix + "charge"         + Suffix );
  produces <std::vector<double> >         ( Prefix + "path"           + Suffix );
  produces <std::vector<double> >         ( Prefix + "thickness"      + Suffix );
  produces <std::vector<double> >         ( Prefix + "localpitch"     + Suffix );
  #ifdef ExtendedCALIBTree
  produces <std::vector<double> >         ( Prefix + "chargeoverpath" + Suffix );
  #endif
  produces <std::vector<unsigned char> >  ( Prefix + "amplitude"      + Suffix );
  produces <std::vector<double> >         ( Prefix + "gainused"       + Suffix );
  produces <std::vector<double> >         ( Prefix + "gainusedTick"   + Suffix );
}

void ShallowGainCalibration::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto         trackindex    = std::make_unique<std::vector<int>>            ();
  auto         rawid         = std::make_unique<std::vector<unsigned int>>   ();
  auto         subdetid      = std::make_unique<std::vector<unsigned int>>   ();
  auto         clusteridx    = std::make_unique<std::vector<unsigned int>>   ();
  auto         stripidx      = std::make_unique<std::vector<unsigned int>>   ();
  auto         localphi      = std::make_unique<std::vector<double>>         ();
  auto         localtheta    = std::make_unique<std::vector<double>>         ();
  auto         localdirx     = std::make_unique<std::vector<double>>         ();
  auto         localdiry     = std::make_unique<std::vector<double>>         ();
  auto         localdirz     = std::make_unique<std::vector<double>>         ();
  auto         firststrip    = std::make_unique<std::vector<unsigned short>> ();
  auto         nstrips       = std::make_unique<std::vector<unsigned short>> ();
  auto         saturation    = std::make_unique<std::vector<bool>>           ();
  auto         overlapping   = std::make_unique<std::vector<bool>>           ();
  auto         farfromedge   = std::make_unique<std::vector<bool>>           ();
  auto         charge        = std::make_unique<std::vector<unsigned int>>   ();
  auto         path          = std::make_unique<std::vector<double>>         ();
  auto         thick         = std::make_unique<std::vector<double>>         ();
  auto         pitch         = std::make_unique<std::vector<double>>         ();
  auto         wheel         = std::make_unique<std::vector<int>>            ();
  #ifdef ExtendedCALIBTree
  auto         chargeoverpath= std::make_unique<std::vector<double>>         ();
  #endif
  auto         amplitude     = std::make_unique<std::vector<unsigned char>>  ();
  auto         gainused      = std::make_unique<std::vector<double>>         ();
  auto         gainusedTick  = std::make_unique<std::vector<double>>         ();

  edm::ESHandle<TrackerGeometry> theTrackerGeometry;         iSetup.get<TrackerDigiGeometryRecord>().get( theTrackerGeometry );  
  m_tracker=&(* theTrackerGeometry );
  edm::ESHandle<SiStripGain> gainHandle;                     iSetup.get<SiStripGainRcd>().get(gainHandle);
  edm::Handle<edm::View<reco::Track> > tracks;	             iEvent.getByToken(tracks_token_, tracks);	  
  edm::Handle<TrajTrackAssociationCollection> associations;  iEvent.getByToken(association_token_, associations);

  uint32_t Idx = 0;
  for( TrajTrackAssociationCollection::const_iterator association = associations->begin(); association != associations->end(); association++) {
       const Trajectory*  traj  = association->key.get();
       const reco::Track* track = association->val.get();

       vector<TrajectoryMeasurement> measurements = traj->measurements();
       for(vector<TrajectoryMeasurement>::const_iterator measurement_it = measurements.begin(); measurement_it!=measurements.end(); measurement_it++){
          TrajectoryStateOnSurface trajState = measurement_it->updatedState();
          if( !trajState.isValid() ) continue;     

          const TrackingRecHit*         hit                 = (*measurement_it->recHit()).hit();
          const SiStripRecHit1D*        sistripsimple1dhit  = dynamic_cast<const SiStripRecHit1D*>(hit);
          const SiStripRecHit2D*        sistripsimplehit    = dynamic_cast<const SiStripRecHit2D*>(hit);
          const SiStripMatchedRecHit2D* sistripmatchedhit   = dynamic_cast<const SiStripMatchedRecHit2D*>(hit);
          const SiPixelRecHit*          sipixelhit          = dynamic_cast<const SiPixelRecHit*>(hit);

          const SiPixelCluster*   PixelCluster = nullptr;
          const SiStripCluster*   StripCluster = nullptr;
          uint32_t                DetId    = 0;
          Idx += 1;

          for(unsigned int h=0;h<2;h++){
            if(!sistripmatchedhit && h==1){
	       continue;
            }else if(sistripmatchedhit  && h==0){
               StripCluster = &sistripmatchedhit->monoCluster();
	           DetId        = sistripmatchedhit->monoId();
            }else if(sistripmatchedhit  && h==1){
               StripCluster = &sistripmatchedhit->stereoCluster();;
	           DetId        = sistripmatchedhit->stereoId();
            }else if(sistripsimplehit){
               StripCluster = (sistripsimplehit->cluster()).get();
    	       DetId        = sistripsimplehit->geographicalId().rawId();
            }else if(sistripsimple1dhit){
               StripCluster = (sistripsimple1dhit->cluster()).get();
    	       DetId        = sistripsimple1dhit->geographicalId().rawId();
            }else if(sipixelhit){
               PixelCluster = (sipixelhit->cluster()).get();
               DetId    = sipixelhit->geographicalId().rawId();
            }else{
               continue;
            }

            SiStripDetId subdet(DetId);
            uint32_t                SubDetId       = subdet.subDetector();
            if( SubDetId<=2 ) continue;

            LocalVector             trackDirection = trajState.localDirection();
            LocalPoint              hitLocalPos    = trajState.localPosition();

/*
            const uint32_t id = hit->detUnit();
            const StripGeomDetUnit* theStripDet = dynamic_cast<const StripGeomDetUnit*>( theTrackerGeometry->idToDet( id ) );
            const float hitstrip_ = theStripDet->specificTopology().strip( hitLocalPos );

            shallow::CLUSTERMAP clustermap = shallow::make_cluster_map(iEvent, clusters_token_);
            shallow::CLUSTERMAP::const_iterator cluster = match_cluster( id, driftedstrip_, clustermap, *clusters);
            if(cluster != clustermap.end()) {
                unsigned i = cluster->second;
                int strip = hitstrip_;
                double localtheta = hit.thetaAtEntry();
                double localphi = hit.phiAtEntry();
                }
            }
*/

            double                  ltheta         = trackDirection.theta();
            double                  lphi           = trackDirection.phi();
            double                  cosine         = trackDirection.z()/trackDirection.mag();
            bool                    Saturation     = false;
            bool                    Overlapping    = false;
            unsigned int            Charge         = 0;
            double                  Path           = (10.0*thickness(DetId))/fabs(cosine);
            double                  Thick          = thickness(DetId);
            double                  Pitch          = localPitch(DetId, hitLocalPos);
            int                     LayerWheel     = layerwheel(&iSetup, DetId, SubDetId);
            double                  PrevGain       = -1;
            double                  PrevGainTick   = -1;
            int                     FirstStrip     = 0;
            int                     NStrips        = 0;

            if(StripCluster){
               const auto           &  Ampls          = StripCluster->amplitudes();
   	                                   FirstStrip     = StripCluster->firstStrip();
                                       NStrips        = Ampls.size();
               int                     APVId          = FirstStrip/128;
//               int                     theChannel     = StripCluster->channel();
//               int                     StripIdx       = theChannel - FirstStrip;

               if(gainHandle.isValid()){ 
                  PrevGain     =  gainHandle->getApvGain(APVId,gainHandle->getRange(DetId, 1),1); 
                  PrevGainTick =  gainHandle->getApvGain(APVId,gainHandle->getRange(DetId, 0),1);           
               }

               for(unsigned int a=0;a<Ampls.size();a++){               
                  Charge+=Ampls[a];
                  if(Ampls[a] >=254)Saturation =true;
                  amplitude->push_back( Ampls[a] );
                  stripidx->push_back( Idx );
               }

               if(FirstStrip==0                                  )Overlapping=true;
               if(FirstStrip==128                                )Overlapping=true;
               if(FirstStrip==256                                )Overlapping=true;
               if(FirstStrip==384                                )Overlapping=true;
               if(FirstStrip==512                                )Overlapping=true;
               if(FirstStrip==640                                )Overlapping=true;

               if(FirstStrip<=127 && FirstStrip+Ampls.size()>127)Overlapping=true;
               if(FirstStrip<=255 && FirstStrip+Ampls.size()>255)Overlapping=true;
               if(FirstStrip<=383 && FirstStrip+Ampls.size()>383)Overlapping=true;
               if(FirstStrip<=511 && FirstStrip+Ampls.size()>511)Overlapping=true;
               if(FirstStrip<=639 && FirstStrip+Ampls.size()>639)Overlapping=true;

               if(FirstStrip+Ampls.size()==127                   )Overlapping=true;
               if(FirstStrip+Ampls.size()==255                   )Overlapping=true;
               if(FirstStrip+Ampls.size()==383                   )Overlapping=true;
               if(FirstStrip+Ampls.size()==511                   )Overlapping=true;
               if(FirstStrip+Ampls.size()==639                   )Overlapping=true;
               if(FirstStrip+Ampls.size()==767                   )Overlapping=true;
            }else if(PixelCluster){
               const auto&             Ampls          = PixelCluster->pixelADC();
               int                     FirstRow       = PixelCluster->minPixelRow();
               int                     FirstCol       = PixelCluster->minPixelCol();
               FirstStrip                             = ((FirstRow/80)<<3 | (FirstCol/52)) * 128; //Hack to save the APVId
               NStrips                                = 0;
               Saturation                             = false;
               Overlapping                            = false;

               for(unsigned int a=0;a<Ampls.size();a++){
                  Charge+=Ampls[a];
                  if(Ampls[a] >=254)Saturation =true;
               }
            }

// shallow::CLUSTERMAP::const_iterator cluster = clustermap.find( std::make_pair( DetId, FirstStrip )); // Idx = cluster->second;

//  edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
//  iEvent.getByToken(theClustersToken_, clusters);
//  edmNew::DetSetVector<SiStripCluster>::const_iterator itClusters=clusters->begin();
//  for(;itClusters!=clusters->end();++itClusters){
//      uint32_t id = itClusters->id();
//      const moduleVars moduleV(id, tTopo);
//    for(edmNew::DetSet<SiStripCluster>::const_iterator cluster=itClusters->begin(); cluster!=itClusters->end();++cluster){

//            subdetid      ->push_back( moduleV.subdetid );

	    #ifdef ExtendedCALIBTree
	    double ChargeOverPath = (double)Charge / Path ;
	    #endif

            trackindex    ->push_back( shallow::findTrackIndex(tracks, track) );
            rawid         ->push_back( DetId );
            subdetid      ->push_back( SubDetId );
            clusteridx    ->push_back( Idx );
            localphi      ->push_back( lphi );
            localtheta    ->push_back( ltheta );
            localdirx     ->push_back( trackDirection.x() );
            localdiry     ->push_back( trackDirection.y() );
            localdirz     ->push_back( trackDirection.z() );
            firststrip    ->push_back( FirstStrip );
            nstrips       ->push_back( NStrips );
            saturation    ->push_back( Saturation );
            overlapping   ->push_back( Overlapping );
            farfromedge   ->push_back( StripCluster ? IsFarFromBorder(&trajState,DetId, &iSetup) : true );
            charge        ->push_back( Charge );
            path          ->push_back( Path );
            thick         ->push_back( Thick );
            pitch         ->push_back( Pitch );
            wheel         ->push_back( LayerWheel );
            #ifdef ExtendedCALIBTree
            chargeoverpath->push_back( ChargeOverPath );
            #endif
            gainused      ->push_back( PrevGain );
            gainusedTick  ->push_back( PrevGainTick );
          }
       }
  }

  iEvent.put(std::move(trackindex),    Prefix + "trackindex"    + Suffix );
  iEvent.put(std::move(rawid     ),    Prefix + "rawid"         + Suffix );
  iEvent.put(std::move(subdetid  ),    Prefix + "subdetid"      + Suffix );
  iEvent.put(std::move(clusteridx),    Prefix + "clusteridx"    + Suffix );
  iEvent.put(std::move(stripidx),      Prefix + "stripidx"      + Suffix );
  iEvent.put(std::move(localphi  ),    Prefix + "localphi"      + Suffix );
  iEvent.put(std::move(localtheta),    Prefix + "localtheta"    + Suffix );
  iEvent.put(std::move(localdirx ),    Prefix + "localdirx"     + Suffix );
  iEvent.put(std::move(localdiry ),    Prefix + "localdiry"     + Suffix );
  iEvent.put(std::move(localdirz ),    Prefix + "localdirz"     + Suffix );
  iEvent.put(std::move(firststrip),    Prefix + "firststrip"    + Suffix );
  iEvent.put(std::move(nstrips),       Prefix + "nstrips"       + Suffix );
  iEvent.put(std::move(saturation),    Prefix + "saturation"    + Suffix );
  iEvent.put(std::move(overlapping),   Prefix + "overlapping"   + Suffix );
  iEvent.put(std::move(farfromedge),   Prefix + "farfromedge"   + Suffix );
  iEvent.put(std::move(charge),        Prefix + "charge"        + Suffix );
  iEvent.put(std::move(path),          Prefix + "path"          + Suffix );
  iEvent.put(std::move(thick),         Prefix + "thickness"     + Suffix );
  iEvent.put(std::move(pitch),         Prefix + "localpitch"    + Suffix );
  iEvent.put(std::move(wheel),         Prefix + "layerwheel"    + Suffix );
  #ifdef ExtendedCALIBTree
  iEvent.put(std::move(chargeoverpath),Prefix + "chargeoverpath"+ Suffix );
  #endif
  iEvent.put(std::move(amplitude),     Prefix + "amplitude"     + Suffix );
  iEvent.put(std::move(gainused),      Prefix + "gainused"      + Suffix );
  iEvent.put(std::move(gainusedTick),  Prefix + "gainusedTick"  + Suffix );
}

/*
void ShallowGainCalibration::beginJob(const edm::EventSetup& iSetup)
{
   printf("Befin JOB\n");

   edm::ESHandle<TrackerGeometry> tkGeom;
   iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );
   vector<GeomDet*> Det = tkGeom->dets();

   edm::ESHandle<SiStripGain> gainHandle;
   iSetup.get<SiStripGainRcd>().get(gainHandle);
   if(!gainHandle.isValid()){printf("\n#####################\n\nERROR --> gainHandle is not valid\n\n#####################\n\n");exit(0);}

   for(unsigned int i=0;i<Det.size();i++){
      DetId  Detid  = Det[i]->geographicalId();
      int    SubDet = Detid.subdetId();

      if( SubDet == StripSubdetector::TIB ||  SubDet == StripSubdetector::TID ||
          SubDet == StripSubdetector::TOB ||  SubDet == StripSubdetector::TEC  ){

          StripGeomDetUnit* DetUnit     = dynamic_cast<StripGeomDetUnit*> (Det[i]);
          if(!DetUnit)continue;

          const StripTopology& Topo     = DetUnit->specificTopology();
          unsigned int         NAPV     = Topo.nstrips()/128;
   
          for(unsigned int j=0;j<NAPV;j++){
                stAPVGain* APV = new stAPVGain;
                APV->DetId         = Detid.rawId();
                APV->APVId         = j;
                APV->PreviousGain  = 1;

                APVsCollOrdered.push_back(APV);
                APVsColl[(APV->DetId<<3) | APV->APVId] = APV;
          }
      }
   }
}


void ShallowGainCalibration::beginRun(edm::Run &, const edm::EventSetup &iSetup){
    printf("BEFIN RUN\n");

    edm::ESHandle<SiStripGain> gainHandle;
    iSetup.get<SiStripGainRcd>().get(gainHandle);
    if(!gainHandle.isValid()){printf("\n#####################\n\nERROR --> gainHandle is not valid\n\n#####################\n\n");exit(0);}

    for(std::vector<stAPVGain*>::iterator it = APVsCollOrdered.begin();it!=APVsCollOrdered.end();it++){
       stAPVGain* APV = *it;
       SiStripApvGain::Range detGainRange = gainHandle->getRange(APV->DetId);
       APV->PreviousGain = *(detGainRange.first + APV->APVId);
    }
}
*/

bool ShallowGainCalibration::IsFarFromBorder(TrajectoryStateOnSurface* trajState, const uint32_t detid, const edm::EventSetup* iSetup)
{ 
  edm::ESHandle<TrackerGeometry> tkGeom; iSetup->get<TrackerDigiGeometryRecord>().get( tkGeom );

  LocalPoint  HitLocalPos   = trajState->localPosition();
  LocalError  HitLocalError = trajState->localError().positionError() ;

  const GeomDetUnit* it = tkGeom->idToDetUnit(DetId(detid));
  if (dynamic_cast<const StripGeomDetUnit*>(it)==nullptr && dynamic_cast<const PixelGeomDetUnit*>(it)==nullptr) {
     std::cout << "this detID doesn't seem to belong to the Tracker" << std::endl;
     return false;
  }

  const BoundPlane plane = it->surface();
  const TrapezoidalPlaneBounds* trapezoidalBounds( dynamic_cast<const TrapezoidalPlaneBounds*>(&(plane.bounds())));
  const RectangularPlaneBounds* rectangularBounds( dynamic_cast<const RectangularPlaneBounds*>(&(plane.bounds())));

  double DistFromBorder = 1.0;    
  double HalfLength     = it->surface().bounds().length() /2.0;

  if(trapezoidalBounds)
  {
      std::array<const float, 4> const & parameters = (*trapezoidalBounds).parameters();
     HalfLength     = parameters[3];
  }else if(rectangularBounds){
     HalfLength     = it->surface().bounds().length() /2.0;
  }else{return false;}

  if (fabs(HitLocalPos.y())+HitLocalError.yy() >= (HalfLength - DistFromBorder) ) return false;

  return true;
}


double ShallowGainCalibration::thickness(DetId id)
{
 map<DetId,double>::iterator th=m_thicknessMap.find(id);
 if(th!=m_thicknessMap.end())
   return (*th).second;
 else {
   double detThickness=1.;
   //compute thickness normalization
   const GeomDetUnit* it = m_tracker->idToDetUnit(DetId(id));
   bool isPixel = dynamic_cast<const PixelGeomDetUnit*>(it)!=nullptr;
   bool isStrip = dynamic_cast<const StripGeomDetUnit*>(it)!=nullptr;
   if (!isPixel && ! isStrip) {
   //FIXME throw exception
      edm::LogWarning("DeDxHitsProducer") << "\t\t this detID doesn't seem to belong to the Tracker";
      detThickness = 1.;
  }else{
      detThickness = it->surface().bounds().thickness();
  }

   m_thicknessMap[id]=detThickness;//computed value
   return detThickness;
 }

}

double ShallowGainCalibration::localPitch(DetId id, LocalPoint p)
{
 map<DetId,double>::iterator th=m_localPitchMap.find(id);
 if(th!=m_localPitchMap.end())
   return (*th).second;
 else {
   double detlocalPitch=1.;
   const GeomDetUnit* it = m_tracker->idToDetUnit(DetId(id));
   bool isPixel = dynamic_cast<const PixelGeomDetUnit*>(it)!=nullptr;
   bool isStrip = dynamic_cast<const StripGeomDetUnit*>(it)!=nullptr;
   if (!isPixel && ! isStrip) {
   //FIXME throw exception
      edm::LogWarning("DeDxHitsProducer") << "\t\t this detID doesn't seem to belong to the Tracker";
      detlocalPitch = 1.;
  }else{
      const StripTopology& topol = (const StripTopology&)it->topology();
      detlocalPitch = topol.localPitch(p);
  }

   m_localPitchMap[id]=detlocalPitch;//computed value
   return detlocalPitch;
 }

}

int ShallowGainCalibration::layerwheel( const edm::EventSetup* iSetup, uint32_t detid, uint32_t subdetid) {
  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup->get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();
  int wheel;
//  StripGeomDetUnit* DetUnit  = dynamic_cast<StripGeomDetUnit*> (Detid);
//  const StripTopology& tTopo = DetUnit->specificTopology();
  if( SiStripDetId::TIB == subdetid ) {
    wheel    = tTopo->tibLayer(detid);
  } else
  if( SiStripDetId::TID == subdetid ) {
    wheel    = tTopo->tidWheel(detid);
  } else
  if( SiStripDetId::TOB == subdetid ) {
    wheel    = tTopo->tobLayer(detid);
  } else
  if( SiStripDetId::TEC == subdetid ) {
    wheel    = tTopo->tecWheel(detid);
  } else {
    cout<<"detid "<<detid<<" subdetid "<<subdetid<<" TIB "<<SiStripDetId::TIB<<" TID "<<SiStripDetId::TID<<" TOB "<<SiStripDetId::TOB<<" TEC "<<SiStripDetId::TEC<<std::endl;
    wheel=-1;
  }

   return wheel;

}


