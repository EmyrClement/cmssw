#ifndef TPGFEReader_h
#define TPGFEReader_h

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

namespace TPGFEReader{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ECONDReader {
  public:    
    ECONDReader(TPGFEConfiguration::Configuration& cfgs) : configs(cfgs), r(0x0), scanMode(false), inspectEvent(0), nShowEvents(0), isTOTUP(false), eventId(0) { initId();}
    ~ECONDReader() {terminate();}
    
    void setModulePath(uint32_t zs, uint32_t sect, uint32_t lnk, uint32_t SiorSci, uint32_t econt_index, uint32_t LDorHD, uint32_t mod_index){
      zside = zs; sector = sect; link = lnk; det = SiorSci;
      econt = econt_index; selTC4 = LDorHD; module = mod_index;
    }
    void setConfigs(TPGFEConfiguration::Configuration& cfgs) {configs = cfgs;}
    void setTotUp(bool totup) {isTOTUP = totup;}
    
    void initId(){
      zside = 0, sector = 0, link = 0, det = 0;
      econt = 0, selTC4 = 1, module = 0, rocn = 0, half = 0;
    }
    
    void checkEvent(uint64_t event) {inspectEvent = event; scanMode = true;}
    void showFirstEvents(uint32_t events) {nShowEvents = events;}
    bool getCheckMode() {return scanMode;}
    uint64_t getCheckedEvent() {return inspectEvent;}
    bool getTotUp() {return isTOTUP;}
    
    void init(uint32_t, uint32_t, uint32_t);
    void getEvents(std::vector<uint64_t>&, uint64_t&, uint64_t&, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>&, std::vector<uint64_t>&);
    void terminate();
    
  private:
    void getModuleData(uint64_t, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>&, std::vector<uint64_t>&);
    // print all (64-bit) words in event
    void event_dump(const Hgcal10gLinkReceiver::RecordRunning *rEvent){
      const uint64_t *p64(((const uint64_t*)rEvent)+1);
      for(uint32_t i(0);i<rEvent->payloadLength();i++){
	std::cout << "EventId : " << std::setw(3) << rEvent->slinkBoe()->eventId() ;
	std::cout << "\t Word " << std::setw(3) << i << " ";
	std::cout << std::hex << std::setfill('0');
	std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	std::cout << std::dec << std::setfill(' ');
      }
      std::cout << std::endl;
    }
    
    bool is_ffsep_word(const uint64_t word, int& isMSB) {

      if( (word & 0xFFFFFFFF) == 0xFFFFFFFF){
	isMSB = 0;
	return true; //0xFFFF 0xFFFF
      }else if((word >> 32) == 0xFFFFFFFF){
	isMSB = 1;
	return true; //0xFFFF 0xFFFF
      }else{
	isMSB = -1;
	return false;
      }
    }

    // returns location in this event of the n'th 0xfecafe... line
    int find_ffsep_word(const Hgcal10gLinkReceiver::RecordRunning *rEvent, int& isMSB, int ffsep_word_loc = -1) {
      const uint64_t *p64(((const uint64_t*)rEvent)+1);
      if(scanMode and eventId==inspectEvent)
	std::cout << "\t\t ffsep_word_loc: "<< ffsep_word_loc
		  << ", word: 0x"<< std::hex << std::setfill('0') << std::setw(16)
		  <<p64[ffsep_word_loc]
		  <<std::dec << std::setfill(' ') << std::setw(4)
		  << ", isMSB : "<<isMSB << std::endl;
      if (ffsep_word_loc > 0) {
	if (is_ffsep_word(p64[ffsep_word_loc], isMSB)) {
	  return ffsep_word_loc;
	}else
	  return ffsep_word_loc;
      } else {
	if(scanMode and eventId==inspectEvent)
	  std::cout << "hgc_roc::ECONDReader::find_ffsep_word(): missed" << std::endl;
      }
 
      return 0;
    }

    bool is_empty_word(const uint64_t word, int& isMSB) {
      if(scanMode and eventId==inspectEvent){
	std::cout << std::hex << std::setfill('0');
	std::cout << "0x" << std::setw(16) << word << std::endl;
	std::cout << std::dec << std::setfill(' ');
      }
      if( (word & 0xFFFFFFFF) == 0x0){
	isMSB = 0;
	return true; //0xFFFF 0xFFFF
      }else if((word >> 32) == 0x0){
	isMSB = 1;
	return true; //0xFECAFE 
      }else{
	isMSB = -1;
	return false;
      }
    }

    bool found_empty_word(const Hgcal10gLinkReceiver::RecordRunning *rEvent, int& isMSB, int ffsep_word_loc = -1) {
      const uint64_t *p64(((const uint64_t*)rEvent)+1);
      if (ffsep_word_loc > 0) {
	if (is_empty_word(p64[ffsep_word_loc], isMSB)) {
	  return true;
	}else
	  return false;
      } 
      else {
	if(scanMode and eventId==inspectEvent) std::cout << "hgc_roc::ECONDReader::found_empty_word(): missed" << std::endl;
	return false;
      }
    }

    bool scanMode;
    uint64_t inspectEvent;
    uint32_t nShowEvents;
    bool isTOTUP;
    
    ////////////////////////////////////////
    //access configs
    ////////////////////////////////////////
    TPGFEConfiguration::Configuration& configs;
    ////////////////////////////////////////

    ////////////////////////////////////////
    //file reader
    ////////////////////////////////////////
    Hgcal10gLinkReceiver::FileReader _fileReader;
    Hgcal10gLinkReceiver::RecordT<4095> *r;
    uint64_t eventId;
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    //id definition
    ////////////////////////////////////////
    uint32_t zside, sector, link, det;
    uint32_t econt, selTC4, module, rocn, half;
    TPGFEConfiguration::TPGFEIdPacking pck;
    ////////////////////////////////////////

    std::vector<uint64_t> captureheader_lst;
    std::map<uint32_t,std::vector<uint64_t>> econheader_lst; //
    std::vector<uint32_t> captureheaderpos_lst;
    std::map<uint32_t,std::vector<uint32_t>> econheaderpos_lst;
    std::map<uint32_t,std::vector<uint32_t>> zeropos_lst;    
    ////////////////////////////////////////
    
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::init(uint32_t relayNumber, uint32_t runNumber, uint32_t linkNumber){
    //Make the buffer space for the records
    r = new Hgcal10gLinkReceiver::RecordT<4095>;
    _fileReader.setDirectory(std::string("dat/Relay")+std::to_string(relayNumber));
    _fileReader.openRun(runNumber,linkNumber);  
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::getModuleData(uint64_t nEvents, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>& hrocarray, std::vector<uint64_t>& events){
    
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    const Hgcal10gLinkReceiver::SlinkBoe *boe = rEvent->slinkBoe();      
    const Hgcal10gLinkReceiver::SlinkEoe *eoe = rEvent->slinkEoe();	
    const uint64_t *p64(((const uint64_t*)rEvent)+1);
    //const uint32_t *p32(((const uint32_t*)rEvent)+2);

    std::map<uint32_t,TPGFEConfiguration::ConfigEconD>& econDPar = configs.getEconDPara();
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    std::string modName = modNameMap.at(std::make_tuple(det, selTC4, module));
    const std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>&seqToRocpin = (pck.getDetType()==0)?configs.getSiSeqToROCpin():configs.getSciSeqToROCpin();
    
    //int nRx = 6;
    bool isTOTUP = false;
    uint32_t rocn = 0;
    uint32_t half = 0;
    
    for (const auto&  econhpos : econheaderpos_lst){
      uint32_t icapblk = econhpos.first;
      uint32_t link = (icapblk==0)? 0 : 1 ;
      std::vector<uint32_t> econhloc = econhpos.second;
      for(uint32_t iecond = 0 ; iecond < econhloc.size() ; iecond++){
	uint32_t idx = pck.packModId(zside, sector, link, det, iecond, selTC4, module);
	int nRx = int(econDPar[idx].getNeRx());
	
	if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){
	  std::cout<< std::hex ;
	  std::cout<<"For econ " << iecond << ", word is 0x"<< std::setfill('0') << std::setw(16) <<p64[econhloc.at(iecond)]<<std::endl;
	  std::cout<<"For econ " << iecond << ", zero word is 0x"<< std::setfill('0') << std::setw(16) <<p64[zeropos_lst[icapblk].at(iecond)]<<std::endl;
	  std::cout<< std::dec << std::setfill(' ') ;
	  std::cout<<"econd words position for econ " << iecond << ", is "<< std::setfill('0') << std::setw(2) <<econhloc.at(iecond)<<std::endl;
	  std::cout<<"zero words position for econ " << iecond << ", is "<< std::setfill('0') << std::setw(2) << zeropos_lst[icapblk].at(iecond)  <<std::endl;
	  
	  //for(int iw=econhloc.at(iecond)+1 ; iw <= zeropos_lst[icapblk].at(iecond) and iecond == 0 ; iw++ ){
	  //for(int iw=econhloc.at(iecond)+1 ; iw <= zeropos_lst[icapblk].at(iecond) ; iw++ ){	  
	  // std::cout<<"word-64 "<< iw << std::hex <<"-th : 0x"<< std::setfill('0') << std::setw(16) << p64[iw] ;
	  // std::cout<<",       word-32 0th : 0x"<< std::setfill('0') << std::setw(8) << p32[2*iw+1] ;
	  // std::cout<<", 0x" << std::setfill('0') << std::setw(8) << p32[2*iw] << std::std::endl;
	  // std::cout<< std::dec << std::setfill(' ') ;
	}
	
	int isMSB[12];
	int ffsep_word_loc[12];
	bool is_ff_sep_notfound[12];
	bool has_all_sep_found = true;;
	int ffword_loc = econhloc.at(iecond) + 1;
	int ffsep_length = 0;
	for(int iloc=0;iloc<nRx;iloc++){
	  isMSB[iloc] = (iloc%2==0) ? 0 : 1 ;
	  ffsep_length = (iloc%2==0) ? 20 : 19 ;
	  is_ff_sep_notfound[iloc] = false;
	  if(scanMode and eventId==inspectEvent)
	    std::cout << "\t iloc : " << iloc <<", isMSB[iloc] : " << isMSB[iloc] <<", ffsep_length : " << ffsep_length << ", ffword_loc: " << ffword_loc << std::endl;
	  ffsep_word_loc[iloc] = find_ffsep_word(rEvent, isMSB[iloc], ffword_loc);
	  ffword_loc += ffsep_length ;
	  if(isMSB[iloc] == -1) {
	    is_ff_sep_notfound[iloc] = true;
	    has_all_sep_found = false ;
	  }
	}
	
	int expctd_empty_word_loc = zeropos_lst[icapblk].at(iecond);
	int empty_isMSB = -1;
	bool hasfound_emptry_word = found_empty_word(rEvent, empty_isMSB, expctd_empty_word_loc) ;      
	if(!hasfound_emptry_word) continue;
	if(!has_all_sep_found) continue;
      
	int ichip = 0;
	  
	for(int iloc=0;iloc<nRx;iloc++){ //For full LD/HD modules : 6/12 respectively
	  uint32_t max_sep_word = (iloc==(nRx-1)) ? expctd_empty_word_loc : ffsep_word_loc[iloc+1] ;
	  bool ishalf = (iloc%2==0) ? false : true ;
	  if(isMSB[iloc]==1) max_sep_word++;
	  int index = 0;
	  uint32_t ch = 0;
	  uint32_t rocpin = 0;
	  int adcmL = 0, adcL = 0, toaL = 0, totL = 0;
	  int adcmM = 0, adcM = 0, toaM = 0, totM = 0;
	  rocn = uint32_t(ichip);
	  half = uint32_t(ishalf);
	  // std::cout<<"iloc : " << iloc << std::endl;
	  TPGFEDataformat::HalfHgcrocChannelData chdata[36];
	  for(uint32_t i = ffsep_word_loc[iloc]+1 ; i< max_sep_word ; i++){	  //For 37 seqs for full half ROC
	    
	    uint32_t wordL = 0;
	    uint32_t wordM = 0;
	    uint32_t seqL = 0;		
	    uint32_t seqM = 0;
	    if(isMSB[iloc]==0){
	      wordL = (p64[i] & 0xFFFFFFFF) ;
	      wordM = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      seqM = 2*index;
	      seqL = 2*index+1;		
	    }else if(isMSB[iloc]==1){
	      wordL = (p64[i-1] & 0xFFFFFFFF) ;	  
	      wordM = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      seqL = 2*index;
	      seqM = 2*index+1;		
	    }
	    adcmL = 0; adcL = 0; toaL = 0; totL = 0;
	    adcmM = 0; adcM = 0; toaM = 0; totM = 0;
	      
	    const uint16_t trigflagL = (wordL>>30) & 0x3;
	    if(trigflagL<=1){ //0 or 1
	      adcmL = (wordL>>20) & 0x3FF;
	      adcL = (wordL>>10) & 0x3FF;
	      toaL = wordL & 0x3FF;
	    }else if(trigflagL==2){
	      adcmL = (wordL>>20) & 0x3FF;
	      totL = (wordL>>10) & 0x3FF;
	      toaL = wordL & 0x3FF;	    
	    }else if(trigflagL==3){
	      adcmL = (wordL>>20) & 0x3FF;
	      totL = (wordL>>10) & 0x3FF;
	      toaL = wordL & 0x3FF;
	    }
	    if(totL>>0x9==1)  totL = (totL & 0x1ff) << 0x3 ; //10-bit to 12-bit conversion
	    if(isTOTUP) totL += 0x7;
	      
	    const uint16_t trigflagM = (wordM>>30) & 0x3;
	    if(trigflagM<=1){ //0 or 1
	      adcmM = (wordM>>20) & 0x3FF;
	      adcM = (wordM>>10) & 0x3FF;
	      toaM = wordM & 0x3FF;
	    }else if(trigflagM==2){
	      adcmM = (wordM>>20) & 0x3FF;
	      totM = (wordM>>10) & 0x3FF;
	      toaM = wordM & 0x3FF;	    
	    }else if(trigflagM==3){
	      adcmM = (wordM>>20) & 0x3FF;
	      totM = (wordM>>10) & 0x3FF;
	      toaM = wordM & 0x3FF;
	    }
	    if(totM>>0x9==1)  totM = (totM & 0x1ff) << 0x3 ; //10-bit to 12-bit conversion
	    if(isTOTUP) totM += 0x7;
	      
	    if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){ 
	      if(isMSB[iloc]==0){
		std::cout<<"\ti : "<<i<<", index : "<<index<<std::endl;
		std::cout << std::hex << std::setfill('0');
		std::cout<<"\tWordM : 0x" << std::setw(8) << wordM <<", wordL : 0x" << std::setw(8) << wordL<< std::endl;
		std::cout << std::dec << std::setfill(' ');
		std::cout<<"\tM:(flag,adc,tot,toa) : (" << trigflagM <<", " << adcM << ", " << totM <<", " << toaM << "), \t"
			 <<"L:(flag,adc,tot,toa) : (" << trigflagL <<", " << adcL << ", " << totL <<", " << toaL << ")" << std::endl;
	      }else{
		std::cout<<"\ti : "<<i<<", index : "<<index<<std::endl;
		std::cout << std::hex << std::setfill('0');
		std::cout<<"\tWordL[i-1] : 0x" << std::setw(8) << wordL<<", wordM : 0x" << std::setw(8) << wordM << std::endl;
		std::cout << std::dec << std::setfill(' ');
		std::cout<<"\tL:(flag,adc,tot,toa)[i-1] : (" << trigflagL <<", " << adcL << ", " << totL <<", " << toaL << "), \t"
			 <<"M:(flag,adc,tot,toa) : (" << trigflagM <<", " << adcM << ", " << totM <<", " << toaM << ")" << std::endl;
	      }
	    }
	      
	    if(seqToRocpin.find(std::make_pair(modName,std::make_tuple(rocn,half,seqM)))!=seqToRocpin.end()){
	      uint32_t rocpin = seqToRocpin.at(std::make_pair(modName,std::make_tuple(rocn,half,seqM))) ;
	      ch = rocpin%36 ; //36 for halfroc
	      if(trigflagM>=0x2){
		chdata[ch].setTot(uint16_t(totM),uint16_t(trigflagM));
		//chdata[ch].setAdc(uint16_t(adcmM),0);
	      }else if(trigflagM==0)
		chdata[ch].setAdc(uint16_t(adcM),uint16_t(trigflagM));
	      else if(trigflagM==0x1){
		chdata[ch].setAdc(0,uint16_t(trigflagM)); //uint16_t(adcM)
		//chdata[ch].setAdc(uint16_t(adcM),uint16_t(trigflagM)); //
		//chdata[ch].setTot(uint16_t(toaM),uint16_t(trigflagM));
	      }else
		chdata[ch].setZero();
	    }
	    if(seqToRocpin.find(std::make_pair(modName,std::make_tuple(rocn,half,seqL)))!=seqToRocpin.end()){
	      uint32_t rocpin = seqToRocpin.at(std::make_pair(modName,std::make_tuple(rocn,half,seqL))) ;
	      ch = rocpin%36 ; //36 for halfroc
	      if(trigflagL>=0x2){
		chdata[ch].setTot(uint16_t(totL),uint16_t(trigflagL));
		//chdata[ch].setAdc(uint16_t(adcmL),0);
	      }else if(trigflagL==0)
		chdata[ch].setAdc(uint16_t(adcL),uint16_t(trigflagL));
	      else if(trigflagL==0x1){
		chdata[ch].setAdc(0,uint16_t(trigflagL)); //uint16_t(adcL)
		//chdata[ch].setAdc(uint16_t(adcL),uint16_t(trigflagL)); //
		//chdata[ch].setTot(uint16_t(toaL),uint16_t(trigflagL));
	      }else
		chdata[ch].setZero();
	    }
	    index++;
	  }//end of loop for 37 seq channels
	  TPGFEDataformat::HalfHgcrocData hrocdata;
	  hrocdata.setChannels(chdata);
	  // uint16_t bxdiff = 3;
	  // int bx = uint16_t((econheader_lst[icapblk].at(iecond)>>20) & 0xFFF) - bxdiff;
	  // if(bx<0) bx = (3564+bx) + 1;
	  uint16_t bx = uint16_t((econheader_lst[icapblk].at(iecond)>>20) & 0xFFF) ;
	  hrocdata.setBx(bx);
	  hrocdata.setSlinkBx(eoe->bxId());
	  pck.setZero();
	  hrocarray[boe->eventId()].push_back(std::make_pair(pck.packRocId(zside, sector, link, det, iecond, selTC4, module, rocn, half),hrocdata));
	  //std::cout<<std::endl;
	  if(iloc%2==1)ichip++;
	}//loop over 6/12 eRxs
	//}//iw loop
      }//econd loop
    }//link/captureblock loop
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::getEvents(std::vector<uint64_t>& refEvents, uint64_t& minEventDAQ, uint64_t& maxEventDAQ, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>& hrocarray, std::vector<uint64_t>& events){
    
    //Set up specific records to interpet the formats
    const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
    const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    const uint64_t *p64(((const uint64_t*)rEvent)+1);
    uint64_t nEvents = 0;
    const int maxEcons = 12;
    uint16_t l1atype = 0;
    const uint8_t valid_daq_be_pattern = 0xfe;
    
    const uint64_t maxEvent = (refEvents.size()==0) ? maxEventDAQ : refEvents[refEvents.size()-1] ;
    //Use the fileReader to read the records
    while(_fileReader.read(r) and nEvents<=maxEvent ) {
      //Check the state of the record and print the record accordingly
      if( r->state()==Hgcal10gLinkReceiver::FsmState::Starting){
	if(!(rStart->valid())){
	  std::cerr << " FsmState::Starting validadity fails : rStart->valid() " << rStart->valid() << std::endl;
	  continue;
	}
      }
      
      else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping){
	if(!(rStop->valid())){
	  std::cerr << " FsmState::Stopping validadity fails : rStop->valid() " << rStop->valid() << std::endl;
	  continue;
	}
      }
      //Else we have an event record 
      else{
	
	const Hgcal10gLinkReceiver::SlinkBoe *boe = rEvent->slinkBoe();      
	const Hgcal10gLinkReceiver::SlinkEoe *eoe = rEvent->slinkEoe();
	const Hgcal10gLinkReceiver::BePacketHeader *beh = rEvent->bePacketHeader();
	
	
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){ 
	  rEvent->RecordHeader::print();
	  event_dump(rEvent);
	  boe->print();
	  eoe->print();
	  std::cout<<"Payload length : "<< rEvent->payloadLength() << ", DAQ-data bxId " << eoe->bxId()  << std::endl;
	}
	
	//Increment event counter
	nEvents++;
	
	if(boe->boeHeader()!=boe->BoePattern) {
	  std::cerr <<"Event: " << nEvents<< " Slink BoE header mismatch " << std::endl;
	  continue;
	}
	if(eoe->eoeHeader()!=eoe->EoePattern) {
	  std::cerr <<"Event: " << nEvents<< " Slink EoE header mismatch " << std::endl;
	  continue;
	}
	
	// if(boe->eventId()%10000==0){
	//   std::cout <<"DAQ: Event: " << boe->eventId() << std::endl;
	//   boe->print();
	//   beh->print();
	//   std::cout <<"=========================================================="<< std::endl;
	// }
	
	if( beh->pattern() != valid_daq_be_pattern ){
	  std::cerr <<"Event: " << boe->eventId() << " DAQ BePacketHeader pattern is  " << std::hex << std::setfill('0') << std::setw(2) << uint32_t(beh->pattern()) << std::dec << std::endl;
	  beh->print();
	  continue;
	}
	
	eventId = boe->eventId();

	l1atype = boe->l1aType();
	
	// if(l1atype==0) {
	//   std::cerr << " Event: "<< eventId << ", l1aType : " << l1atype << std::endl;
	//   //Event_Dump(eventId, rEvent);
	//   continue;
	// }

	bool eventscan = (refEvents.size()==0) ? (boe->eventId()>=minEventDAQ and boe->eventId()<maxEventDAQ) : (std::find(refEvents.begin(),refEvents.end(),eventId) != refEvents.end()) ;
	// /////////////////////////////////////////////////////////////////
	//if(boe->eventId()>=minEventDAQ and boe->eventId()<maxEventDAQ){
	if(eventscan){
	  captureheader_lst.clear();
	  econheader_lst.clear();
	  captureheaderpos_lst.clear();
	  econheaderpos_lst.clear();
	  zeropos_lst.clear();
	  
	  uint32_t icapblk = 0;
	  uint32_t capturepos = 2;
	  bool isThreeEcons = true;
	  bool isallPassthrough = true;
	  bool isValideRxSize = true;
	  while( (capturepos+2) < rEvent->payloadLength() ){
	    uint64_t captureblkheader = p64[capturepos];
	    captureheader_lst.push_back(captureblkheader);
	    captureheaderpos_lst.push_back(capturepos);
      
	    bool isEcondp[maxEcons];
	    bool isPassthrough[maxEcons];
	    uint32_t zeroloc[maxEcons];
	    int nofecons = 0;
	    for(int iecon=0;iecon<maxEcons;iecon++) {
	      isEcondp[iecon] = false;
	      isPassthrough[iecon]  = false;
	      zeroloc[iecon] = 0;
	    }
	    for(int iecon=0;iecon<maxEcons;iecon++){
	      uint8_t stat = (p64[capturepos] >> uint8_t(iecon*3)) & 0x7;
	      if(stat==0b000 or stat==0b010 or stat==0b011) {
		isEcondp[iecon] = true ;
		nofecons++;
	      }
	    }
	    if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)) std::cout << "icapblk: " << icapblk << ", nofecons: " << nofecons << std::endl;
	    //assert(nofecons!=0);
	    if(nofecons!=3) {
	      isThreeEcons = false;
	      break;
	     } //update the feedback from Martim
	    uint32_t next_econd_pos = capturepos + 1;
	    for(int iecond = 0 ; iecond < nofecons ; iecond++){
	      uint32_t econpos = next_econd_pos;
	      uint32_t econh0 = (p64[econpos] >> 32) & 0xFFFFFFFF ;
	      uint32_t econh1 = p64[econpos] & 0xFFFFFFFF ;
	      econheader_lst[icapblk].push_back(p64[econpos]); 
	      econheaderpos_lst[icapblk].push_back(econpos);
	      isPassthrough[iecond] = (econh0 >> 13 ) & 0x1 ;
	      //assert(isPassthrough[iecond]);
	      if(!isPassthrough[iecond]) {
		isallPassthrough = false;
		break;
	      }
	      uint32_t econsize = (econh0 >> 14 ) & 0x1FF ; //specifies length in 32b-words+1 of eRx sub-packets and CRC trailer
	      //assert(econsize>1);
	      if((econsize-1)<=0){
		isValideRxSize = false;
		break;
	      }
	      zeroloc[iecond] = (econsize-1)/2 + (econpos+1) ; //(econ0pos+1): +1 to start from first eRx header
	      zeropos_lst[icapblk].push_back(zeroloc[iecond]);
	      next_econd_pos = zeroloc[iecond]+1 ;
	      if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){
		std::cout << "icapblk: " << icapblk << ", iecond: " << iecond << ", econh0: " << econh0 << ", econh1: " << econh1
			  <<", econsize: " << econsize
			  << ", zeroloc[iecond]: " << zeroloc[iecond] << std::endl;
	      }
	    }
	    capturepos = zeroloc[nofecons-1] + 1 ;
	    icapblk++;
	  }//find the block positions
	  if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){
	    std::cout<< " icapblk: " << icapblk << std::endl;
	    std::cout<< " econheaderpos_lst.size(): " << econheaderpos_lst.size() << std::endl;
	    std::cout<< " isThreeEcons: "<< isThreeEcons << ", isallPassthrough: "<< isallPassthrough << ", isValideRxSize: " << isValideRxSize << std::endl;
	    for (const auto&  econhpos : econheaderpos_lst)
	      for(int iecond = 0 ; iecond < econheaderpos_lst[econhpos.first].size() ; iecond++)
		std::cout<<"icapblk: "<<econhpos.first<<" zero position for econ " << iecond << ", "<<zeropos_lst[econhpos.first].at(iecond)<<std::endl;
	  }
	  //if(isThreeEcons and isallPassthrough and isValideRxSize){
	  if(isallPassthrough and isValideRxSize){
	    getModuleData(nEvents, hrocarray, events);
	    events.push_back(boe->eventId());
	  }
	}//read event if within the range
      }//else for loopting events
    }//file reader
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::terminate(){
    _fileReader.close();
    if(!r) delete r;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ECONTReader {
  public:
    
    ECONTReader(TPGFEConfiguration::Configuration& cfgs) : configs(cfgs), r(0x0), scanMode(false), inspectEvent(0), nShowEvents(0), eventId(0), moduleId(0), maxCAFESeps(7) { initId();}
    ~ECONTReader() {terminate();}
    
    void setModulePath(uint32_t zs, uint32_t sect, uint32_t lnk, uint32_t SiorSci, uint32_t econt_index, uint32_t LDorHD, uint32_t mod_index){
      zside = zs; sector = sect; link = lnk; det = SiorSci;
      econt = econt_index; selTC4 = LDorHD; module = mod_index;
    }
    void setConfigs(TPGFEConfiguration::Configuration& cfgs) {configs = cfgs;}
    void initId(){
      zside = 0, sector = 0, link = 0, det = 0;
      econt = 0, selTC4 = 1, module = 0, rocn = 0, half = 0;
    }

    void checkEvent(uint64_t event) {inspectEvent = event; scanMode = true;}
    void showFirstEvents(uint32_t events) {nShowEvents = events;}
    bool getCheckMode() {return scanMode;}
    uint64_t getCheckedEvent() {return inspectEvent;}
    
    void setNofCAFESep(int maxsep) { maxCAFESeps = maxsep;}
    void init(uint32_t, uint32_t, uint32_t);
        
    void getEvents(std::vector<uint64_t>& refEvents, uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGBEDataformat::Trig24Data>>>& econtarray,
		   std::map<uint64_t,std::vector<std::pair<uint32_t,TPGBEDataformat::TrigTCProcData>>>& tcprocarray);
    void terminate();
    
  private:
    // print all (64-bit) words in event
    void event_dump(const Hgcal10gLinkReceiver::RecordRunning *rEvent){
      const uint64_t *p64(((const uint64_t*)rEvent)+1);
      for(uint32_t i(0);i<rEvent->payloadLength();i++){
	std::cout << "EventId : " << std::setw(3) << rEvent->slinkBoe()->eventId() ;
	std::cout << "\t Word " << std::setw(3) << i << " ";
	std::cout << std::hex << std::setfill('0');
	std::cout << "0x" << std::setw(16) << p64[i] << std::endl;
	std::cout << std::dec << std::setfill(' ');
      }
      std::cout << std::endl;
    }
    
    uint64_t Abs32(uint32_t a, uint32_t b){
      if(a>b)
	return (a-b);
      else if(b>a)
	return (b-a);
      else
	return 0;
    }
    
    uint64_t Abs64(uint64_t a, uint64_t b){
      if(a>b)
	return (a-b);
      else if(b>a)
	return (b-a);
      else
	return 0;
    }

    bool is_cafe_word(const uint64_t word) {
      if(((word >> 32) & 0xFFFFFFFF) == 0xcafecafe)
	return true; //0xFECAFE                                                      // 0xFF FFFF masking to be there to compare with 16698110 (could have used the hex of this)
      else
	return false;
    }   

    // returns location in this event of the n'th 0xfecafe... line
    int find_cafe_word(const Hgcal10gLinkReceiver::RecordRunning *rEvent, int n, int cafe_word_loc = -1) {
      const uint64_t *p64(((const uint64_t*)rEvent)+1);

      if (cafe_word_loc > 0) {
	if (is_cafe_word(p64[cafe_word_loc])) {
	  return cafe_word_loc;
	}else
	  return cafe_word_loc;
      } 
      else {
	;;//std::cout << "oops" << std::endl;
      }
 
      int cafe_counter = 0;
      int cafe_word_idx = -1;
      for(unsigned i(0);i<rEvent->payloadLength();i++){
	const uint64_t word = p64[i];                                                                                // Is it advisable to convert 64 bit to 32 this way without masking ?
	if (is_cafe_word(word)) { // if word == 0xfeca
	  cafe_counter++;       
	  if (cafe_counter == n){                                                                                    
	    cafe_word_idx = i;
	    break;
	  }
	}
      }

      if (cafe_word_idx == -1) {
	//std::cerr << "Could not find cafe word" << std::endl;
	return 0;
      }else {
	return cafe_word_idx;
      }
    }
    
    void set_packet_tcs_bc(bool packet_tcs[48], uint32_t* packet, const uint32_t& maxNofTcs) {
      //The first 19 bits are in the first word
      // for (int i=0; i<=19; i++) {
      //   int j = 47 - i;
      //   packet_tcs[j] = (pick_bits32(packet[0], 12+i, 1)==1)?true:false;
      // }
      // for (int i=0; i<=27; i++) {
      //   int j = 27 - i;
      //   packet_tcs[j] = (pick_bits32(packet[1], i, 1)==1)?true:false;
      // }
      
      //Will it be always bit-pack for 48 Tcs or less for partial modules
      for (int i=0; i<=19; i++) {
	packet_tcs[i] = (pick_bits32(packet[0], 12+i, 1)==1)?true:false;
      }
      for (int i=0; i<=27; i++) {
	int j = 20 + i;
	packet_tcs[j] = (pick_bits32(packet[1], i, 1)==1)?true:false;
      }
    }

    // 9 energies, 7 bits long, immediately after the packet energies
    void set_packet_energies_bc9(uint64_t packet_energies[9], uint32_t* packet, const uint32_t& maxBCTcs) {
      
      uint64_t packet64[4];
      for (int i=0; i<4; i++) {
	packet64[i] = packet[i];
      }

      // based on the number of BC Tcs the following can be rearranged. In addition, the size of above "packet" array might change as well.
      
      // need one 64 bit words since all of the energies are 9*7 = 63 bits long
      uint64_t word1 = (packet64[1] << (28+32)) + (packet64[2] << 28) + (packet64[3] >> 4);
  
      for (int i=0; i<9; i++) 
	packet_energies[i] = pick_bits64(word1, i*7, 7);
    }

    uint32_t pick_bits32(uint32_t number, int start_bit, int number_of_bits) {
      // Create a mask to extract the desired bits.
      uint32_t mask = (1 << number_of_bits) - 1;
      // Shift the number to the start bit position.
      number = number >> (32 - start_bit - number_of_bits);
      // Mask the number to extract the desired bits.
      uint32_t picked_bits = number & mask;

      return picked_bits;
    }
    
    uint64_t pick_bits64(uint64_t number, int start_bit, int number_of_bits) {
      // Create a mask to extract the desired bits.
      uint64_t mask = (1 << number_of_bits) - 1;
      // Shift the number to the start bit position.
      number = number >> (64 - start_bit - number_of_bits);
      // Mask the number to extract the desired bits.
      uint64_t picked_bits = number & mask;

      return picked_bits;
    }

    // 12 locations, 2 bits long, immediately after the packet counter
    void set_packet_locations(uint32_t packet_locations[12], uint32_t* packet) {
      for (int i=0; i<12; i++) {
	packet_locations[i] = pick_bits32(packet[0], 4+i*2, 2);
      }
    }
    
    // 12 energies, 7 bits long, immediately after the packet energies
    void set_packet_energies(uint64_t packet_energies[12], uint32_t* packet) {
      uint64_t packet64[4];
      for (int i=0; i<4; i++) {
	packet64[i] = packet[i];
      }

      // need two 64 bit words since all of the energies are 12*7 = 84 bits long
      // word 1 starts with the beginning of the energies
      uint64_t word1 = (packet64[0] << (28+32)) + (packet64[1] << 28) + (packet64[2] >> 4);
      // word 2 are the last 64 bits of the packet (which is 128 bits long)
      uint64_t word2 = (packet64[2] << 32) + packet64[3];

      for (int i=0; i<12; i++) {
	if (i < 9) {
	  // first 9 (0->8) energies fit in first word
	  packet_energies[i] = pick_bits64(word1, i*7, 7);
	} 
	else {
	  // 9th energy starts 27 bits into the second word
	  packet_energies[i] = pick_bits64(word2, 27+(i-9)*7, 7); 
	} 
      }
    }
    
    ////////////////////////////////////////
    bool scanMode;
    uint64_t inspectEvent;
    uint32_t nShowEvents;
    bool isMSB;
    int maxCAFESeps;
    ////////////////////////////////////////
    //access configs
    ////////////////////////////////////////
    TPGFEConfiguration::Configuration& configs;
    ////////////////////////////////////////

    ////////////////////////////////////////
    //file reader
    ////////////////////////////////////////
    Hgcal10gLinkReceiver::FileReader _fileReader;
    Hgcal10gLinkReceiver::RecordT<4095> *r;
    uint64_t eventId;
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    //id definition
    ////////////////////////////////////////
    uint32_t zside, sector, link, det;
    uint32_t econt, selTC4, module, rocn, half;
    TPGFEConfiguration::TPGFEIdPacking pck;
    uint32_t moduleId;
    ////////////////////////////////////////
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONTReader::init(uint32_t relayNumber, uint32_t runNumber, uint32_t linkNumber){
    //Make the buffer space for the records
    // uint32_t trig_linkNumber = TMath::FloorNint((linkNumber-1)/2);
    // r = new Hgcal10gLinkReceiver::RecordT<4095>;
    // _fileReader.setDirectory(std::string("dat/Relay")+std::to_string(relayNumber));
    // _fileReader.openRun(runNumber,trig_linkNumber);
    // isMSB = (linkNumber%2==0) ? true : false;
    
    r = new Hgcal10gLinkReceiver::RecordT<4095>;
    _fileReader.setDirectory(std::string("dat/Relay")+std::to_string(relayNumber));
    _fileReader.openRun(runNumber,linkNumber);  

  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  void ECONTReader::getEvents(std::vector<uint64_t>& refEvents, uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGBEDataformat::Trig24Data>>>& econtarray,
			      std::map<uint64_t,std::vector<std::pair<uint32_t,TPGBEDataformat::TrigTCProcData>>>& tcprocarray){

    //Set up specific records to interpet the formats
    const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
    const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    uint64_t nEvents = 0;
    
    uint64_t prevEvent = 0;
    uint32_t prevSequence = 0;
    uint64_t eventId = 0;
    uint32_t sequenceId = 0;
    uint16_t bxId = 0;
    uint16_t l1atype = 0;
    int ievent = 0;

    uint64_t maxEvent = (refEvents.size()==0) ? maxEventTrig : refEvents[refEvents.size()-1] ;
    //Use the fileReader to read the records
    while(_fileReader.read(r) and nEvents<=maxEvent ) {
      //Check the state of the record and print the record accordingly
      if( r->state()==Hgcal10gLinkReceiver::FsmState::Starting){
	if(!(rStart->valid())){
	  std::cerr << " FsmState::Starting validadity fails : rStart->valid() " << rStart->valid() << std::endl;
	  continue;
	}
	
      } else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping){	
	if(!(rStop->valid())){
	  std::cerr << " FsmState::Stopping validadity fails : rStop->valid() " << rStop->valid() << std::endl;
	  continue;
	}
      }
      //Else we have an event record 
      else{
	
	const Hgcal10gLinkReceiver::SlinkBoe *boe = rEvent->slinkBoe();      
	const Hgcal10gLinkReceiver::SlinkEoe *eoe = rEvent->slinkEoe();
	const Hgcal10gLinkReceiver::BePacketHeader *beh =  rEvent->bePacketHeader();
	
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){  
	  rEvent->RecordHeader::print();
	  event_dump(rEvent);
	}
	
	//Increment event counter and reset error state
	nEvents++;
	if(boe->boeHeader()!=boe->BoePattern) {
	  std::cerr <<"Event: " << nEvents<< " Slink BoE header mismatch " << std::endl;
	  continue;
	}
	if(eoe->eoeHeader()!=eoe->EoePattern) {
	  std::cerr <<"Event: " << nEvents<< " Slink EoE header mismatch " << std::endl;
	  continue;
	}

	// if(boe->eventId()%10000==0){
	//   std::cout <<"TRIG: Event: " << boe->eventId() << std::endl;
	//   boe->print();
	//   beh->print();
	//   std::cout <<"=========================================================="<< std::endl;
	// }

	if(uint32_t(beh->pattern())!=0xca){
	  beh->print();
	  std::cerr <<"Event: " << nEvents<< " TPG BePacketHeader pattern is  " << std::hex << std::setfill('0') << std::setw(2) << uint32_t(beh->pattern()) << std::dec << std::endl;
	  continue;
	}
	
	eventId = boe->eventId();
	bxId = eoe->bxId();
	sequenceId = rEvent->RecordHeader::sequenceCounter(); 
	l1atype = boe->l1aType();
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)) std::cout << "L1aType : " << l1atype << ", bxId: " << bxId <<  std::endl;
	if(l1atype==0) {
	  std::cerr << " Event: "<< eventId << ", l1aType : " << l1atype << std::endl;
	  //Event_Dump(eventId, rEvent);
	  continue;
	}
	
	if((Abs64(eventId,prevEvent) != Abs32(sequenceId, prevSequence)) and Abs64(eventId,prevEvent)>=2){
	  std::cerr << " Event: "<< eventId << ", l1aType : " << l1atype << ", and prevEvent  "<< prevEvent << ", nEvents : " << nEvents << " differs by "<< Abs64(eventId,prevEvent) <<" (sequence differs by [ "<< sequenceId << " - "<< prevSequence << " ] = " << Abs32(sequenceId, prevSequence) << "), EventID_Diff is more than 2 " << std::endl;
	  //Event_Dump(eventId, rEvent);
	  rEvent->RecordHeader::print();
	  prevEvent = boe->eventId();
	  prevSequence = rEvent->RecordHeader::sequenceCounter(); 
	  continue;
	  //break;
	}
	
	prevEvent = eventId;
	prevSequence = sequenceId;
	
	int extra_cafe_word_loc = find_cafe_word(rEvent,(maxCAFESeps+1));
	if(extra_cafe_word_loc!=0){
	  std::cerr << " Event: " << eventId << " excess CAFE separator of " <<(maxCAFESeps+1)<< std::endl;
	  continue;
	}
	
	int first_cafe_word_loc = find_cafe_word(rEvent,1);
	if( first_cafe_word_loc != 2){
	  std::cerr << " Event: " << eventId << ", first_cafe_word_loc  "<< first_cafe_word_loc << std::endl;
	  //Event_Dump(eventId, rEvent);
	  continue;
	}
	
	const uint64_t *p64(((const uint64_t*)rEvent)+1);
	
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){ 
	  std::cout<<"========= event : "<< nEvents << "=================="<<  std::endl;
	  std::cout<< std::hex   ;
	  std::cout<<"0th : 0x"<< std::setfill('0') << std::setw(16) << p64[0] <<  std::endl;
	  std::cout<<"1st : 0x" << p64[1] <<  std::endl;
	  std::cout<<"2nd : 0x" << p64[2] <<  std::endl;
	  std::cout<< std::dec << std::setfill(' ') ;
	}
	
	bool eventscan = (refEvents.size()==0) ? (eventId>=minEventTrig and eventId<maxEventTrig) : (std::find(refEvents.begin(),refEvents.end(),eventId) != refEvents.end()) ;
	
	//if(eventId>=minEventTrig and eventId<maxEventTrig){
	if(eventscan){
	  
	  ////////////// First find the block details from header //////////
	  int loc[maxCAFESeps], size[maxCAFESeps];
	  for(int iloc = 0 ; iloc < maxCAFESeps ; iloc++) {loc[iloc] = size[iloc] = -1;}
	  std::cout<< std::dec << std::setfill(' ');
	  for(int iloc = 0 ; iloc < maxCAFESeps ; iloc++){
	    loc[iloc] = find_cafe_word(rEvent, iloc+1);
	    size[iloc] = p64[loc[iloc]] & 0xFF ;
	    int chid = (p64[loc[iloc]] >> 8) & 0xFF ;
	    int bufstat = (p64[loc[iloc]] >> 16) & 0xF ;
	    int nofwd_perbx = (p64[loc[iloc]] >> 20) & 0xF ;
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iloc
		       << std::hex
		       <<", word : 0x" << std::setfill('0') << std::setw(16) << p64[loc[iloc]]
		       << std::dec << std::setfill(' ')
		       << ", location : " << loc[iloc] << ", size: " << size[iloc] <<", ch_id : "<< chid << ", bufstat: " << bufstat << ", nofwd_perbx: "<< nofwd_perbx <<  std::endl;
	  }
	  /////////////////////////////////////////////////////////////////
	  
	  bool isExpectedBlockSize = true;
	  for(int iloc = 0 ; iloc < maxCAFESeps ; iloc++) {
	    bool checksize = false;
	    if(iloc==(maxCAFESeps-1)){
	      int withzeropadding = (maxCAFESeps==14)?4:3;
	      int totsize = loc[iloc] + size[iloc] + withzeropadding ; //zero padding is added to make the event size even
	      checksize = (totsize==rEvent->payloadLength());
	      if(!checksize){
		std::cerr << " Event: "<< eventId << " has mismatch in last cafe position : " << iloc << ", size from block " << totsize << ", payload size " << rEvent->payloadLength() << std::endl;
		// rEvent->RecordHeader::print();
		// boe->print();
		// eoe->print();
		// Event_Dump(eventId, rEvent);
		continue;
	      }
	    }else{
	      int totsize = loc[iloc] + size[iloc] + 1;
	      checksize = (totsize==loc[iloc+1]);
	      if(!checksize){
		std::cerr << " Event: "<< eventId << " has mismatch in cafe position : " << iloc << std::endl;
		//Event_Dump(eventId, rEvent);
		continue;
	      }
	    }
	    if(!checksize) isExpectedBlockSize = false;
	  }
	  if(!isExpectedBlockSize){
	    std::cerr << " Event: "<< eventId << " has block size and location mismatch."<< std::endl;
	    //Event_Dump(eventId, rEvent);
	    continue;
	  }
	  
	  //if(nEvents==1) continue;
	  //if(ievent>(maxEvents-1)) continue;
	  if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)) std::cout<<"iEvent: " << ievent <<  std::endl;
	  
	  //////////// Read raw elink inputs for train 1 /////////////////////
	  int iblock = 1;
	  int elBgnOffset = 0;
	  int elIndx = 0;
	  uint32_t elpckt[5][7][7]; //4:elink-TDAQs + 1 external scintillator TDAQ, the first 7 is for bx and second one for number of elinks
	  uint32_t bx = 0xF;
	  int iel = 0;
	  int ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t wMSB = (p64[iw] >> 32) & 0xFFFFFFFF ;
	    uint32_t wLSB = p64[iw] & 0xFFFFFFFF ;
	    if((elIndx-elBgnOffset)%4==0) iel = 0;
	    if(elIndx>=elBgnOffset){
	      elpckt[0][ibx][iel] = wMSB;
	      if(iel<=6) elpckt[0][ibx][iel+1] = wLSB;
	      iel += 2;
	    }//pick the first elink
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw << ", bx: " << bx << ", ibx: " << ibx
		       << std::hex
		       <<", MSB-word : 0x" << std::setfill('0') << std::setw(8) << wMSB
		       <<", LSB-word : 0x" << std::setfill('0') << std::setw(8) << wLSB
		       << std::dec << std::setfill(' ')
		       << std::endl;
	    elIndx++;
	    if((elIndx-elBgnOffset)%4==0) ibx++;
	  }
	  /////////////////////////////////////////////////////////////////

	  //////////// Read raw elink inputs for train 2 /////////////////////
	  iblock = 2;
	  elBgnOffset = 0;
	  elIndx = 0;
	  bx = 0xF;
	  iel = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t wMSB = (p64[iw] >> 32) & 0xFFFFFFFF ;
	    uint32_t wLSB = p64[iw] & 0xFFFFFFFF ;
	    if((elIndx-elBgnOffset)%4==0) iel = 0;
	    if(elIndx>=elBgnOffset){
	      elpckt[1][ibx][iel] = wMSB;
	      if(iel<=6) elpckt[1][ibx][iel+1] = wLSB;
	      iel += 2;
	    }//pick the first elink
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw << ", bx: " << bx << ", ibx: " << ibx
		       << std::hex
		       <<", MSB-word : 0x" << std::setfill('0') << std::setw(8) << wMSB
		       <<", LSB-word : 0x" << std::setfill('0') << std::setw(8) << wLSB
		       << std::dec << std::setfill(' ')
		       << std::endl;
	    elIndx++;
	    if((elIndx-elBgnOffset)%4==0) ibx++;
	  }
	  /////////////////////////////////////////////////////////////////

	  //////////// Read raw elink inputs for train 3 /////////////////////
	  iblock = 3;
	  elBgnOffset = 0;
	  elIndx = 0;
	  bx = 0xF;
	  iel = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t wMSB = (p64[iw] >> 32) & 0xFFFFFFFF ;
	    uint32_t wLSB = p64[iw] & 0xFFFFFFFF ;
	    if((elIndx-elBgnOffset)%4==0) iel = 0;
	    if(elIndx>=elBgnOffset){
	      elpckt[2][ibx][iel] = wMSB;
	      if(iel<=6) elpckt[2][ibx][iel+1] = wLSB;
	      iel += 2;
	    }//pick the first elink
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw << ", bx: " << bx << ", ibx: " << ibx
		       << std::hex
		       <<", MSB-word : 0x" << std::setfill('0') << std::setw(8) << wMSB
		       <<", LSB-word : 0x" << std::setfill('0') << std::setw(8) << wLSB
		       << std::dec << std::setfill(' ')
		       << std::endl;
	    elIndx++;
	    if((elIndx-elBgnOffset)%4==0) ibx++;
	  }
	  /////////////////////////////////////////////////////////////////
	  
	  //////////// Read raw elink inputs for train 4 /////////////////////
	  iblock = 4;
	  elBgnOffset = 0;
	  elIndx = 0;
	  bx = 0xF;
	  iel = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t wMSB = (p64[iw] >> 32) & 0xFFFFFFFF ;
	    uint32_t wLSB = p64[iw] & 0xFFFFFFFF ;
	    if((elIndx-elBgnOffset)%4==0) iel = 0;
	    if(elIndx>=elBgnOffset){
	      elpckt[3][ibx][iel] = wMSB;
	      if(iel<=6) elpckt[3][ibx][iel+1] = wLSB;
	      iel += 2;
	    }//pick the first elink
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw << ", bx: " << bx << ", ibx: " << ibx
		       << std::hex
		       <<", MSB-word : 0x" << std::setfill('0') << std::setw(8) << wMSB
		       <<", LSB-word : 0x" << std::setfill('0') << std::setw(8) << wLSB
		       << std::dec << std::setfill(' ')
		       << std::endl;
	    elIndx++;
	    if((elIndx-elBgnOffset)%4==0) ibx++;
	  }
	  /////////////////////////////////////////////////////////////////	  
	  
	  //////////// Print unpacker output for train 1 /////////////////////
	  iblock = 5;
	  int unpkBgnOffset = 0;
	  int unpkIndx = 0;
	  uint32_t unpackedWord[5][3][7][2][8]; //4:TDAQ+1 external scintillator TDAQ, 3:ECON-T, 7:bxs, 2:streams, 8:words
	  uint32_t iunpkw = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
	    uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
	    uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
	    uint32_t col3 = p64[iw] & 0xFFFF ;
	    if(unpkIndx>=unpkBgnOffset and (unpkIndx-unpkBgnOffset)%8==0) iunpkw=0;
	    if(unpkIndx>=unpkBgnOffset){
	      unpackedWord[0][0][ibx][0][iunpkw] = col3; //BC6
	      unpackedWord[0][1][ibx][0][iunpkw] = col1; //BC4
	      unpackedWord[0][2][ibx][0][iunpkw] = col0; //BC4
	      iunpkw++;
	    }
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw
		       << std::hex
		       <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
		       <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
		       <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
		       <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
		       << std::dec << std::setfill(' ')
		       << std::endl;
	
	    unpkIndx++;
	    if(unpkIndx>0 and (unpkIndx-unpkBgnOffset)%8==0) ibx++;
	  }
	  /////////////////////////////////////////////////////////////////
	  
	  //////////// Print unpacker output for train 2 /////////////////////
	  iblock = 6;
	  unpkBgnOffset = 0;
	  unpkIndx = 0;
	  iunpkw = 0;
	  ibx = 0;
	  uint32_t col2_prev,col3_prev;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
	    uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
	    uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
	    uint32_t col3 = p64[iw] & 0xFFFF ;
	    if(unpkIndx>=unpkBgnOffset and (unpkIndx-unpkBgnOffset)%8==0) iunpkw=0;
	    if(unpkIndx>=unpkBgnOffset){
	      unpackedWord[1][0][ibx][0][iunpkw] = col3; //STC16
	      unpackedWord[1][1][ibx][0][iunpkw] = col1; //STC16
	      unpackedWord[1][2][ibx][0][iunpkw] = col0; //STC16
	      iunpkw++;
	    }
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw
		       << std::hex
		       <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
		       <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
		       <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
		       <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
		       << std::dec << std::setfill(' ')
		       << std::endl;
		
	    unpkIndx++;
	    if(unpkIndx>0 and (unpkIndx-unpkBgnOffset)%8==0) ibx++;
	  }
	  
	  //////////// Print unpacker output for train 3 /////////////////////
	  iblock = 7;
	  unpkBgnOffset = 0;
	  unpkIndx = 0;
	  iunpkw = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
	    uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
	    uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
	    uint32_t col3 = p64[iw] & 0xFFFF ;
	    if(unpkIndx>=unpkBgnOffset and (unpkIndx-unpkBgnOffset)%8==0) iunpkw=0;
	    if(unpkIndx>=unpkBgnOffset){
	      unpackedWord[2][0][ibx][0][iunpkw] = col3; //STC4-12TCs
	      unpackedWord[2][0][ibx][1][iunpkw] = col2; //STC4-12TCs
	      unpackedWord[2][1][ibx][0][iunpkw] = col1; //STC4-6 //partial module
	      iunpkw++;
	    }
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw
		       << std::hex
		       <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
		       <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
		       <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
		       <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
		       << std::dec << std::setfill(' ')
		       << std::endl;
		
	    unpkIndx++;
	    if(unpkIndx>0 and (unpkIndx-unpkBgnOffset)%8==0) ibx++;
	  }
	  
	  //////////// Print unpacker output for train 4 /////////////////////
	  iblock = 8;
	  unpkBgnOffset = 0;
	  unpkIndx = 0;
	  iunpkw = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
	    uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
	    uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
	    uint32_t col3 = p64[iw] & 0xFFFF ;
	    if(unpkIndx>=unpkBgnOffset and (unpkIndx-unpkBgnOffset)%8==0) iunpkw=0;
	    if(unpkIndx>=unpkBgnOffset){
	      unpackedWord[3][0][ibx][0][iunpkw] = col3; //STC4A-12TCs
	      unpackedWord[3][0][ibx][1][iunpkw] = col2; //STC4A-12TCs
	      unpackedWord[3][1][ibx][0][iunpkw] = col1; //STC4A-10TCs
	      unpackedWord[3][1][ibx][1][iunpkw] = col0; //STC4A-10TCs
	      iunpkw++;
	    }
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
		  std::cout<<"iloc: "<< iw
		       << std::hex
		       <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
		       <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
		       <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
		       <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
		       << std::dec << std::setfill(' ')
		       << std::endl;
		
	    unpkIndx++;
	    if(unpkIndx>0 and (unpkIndx-unpkBgnOffset)%8==0) ibx++;
	  }
	  
	  //////////// Read raw elink like inputs for exterrnal trigger /////////////////////
	  iblock = 9;
	  elIndx = 0;
	  bx = 0xF;
	  iel = 0;
	  ibx = 0;
	  for(int iw = loc[iblock]+1; iw <= (loc[iblock]+size[iblock]) ; iw++ ){
	    uint32_t wLSB = p64[iw] & 0xFFFFFFFF ;
	    elpckt[4][ibx][iel] = wLSB;
	    iel++;	    
	    if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	      std::cout<<"iloc: "<< iw << ", bx: " << bx << ", ibx: " << ibx
		       << std::hex
		       <<", LSB-word : 0x" << std::setfill('0') << std::setw(8) << wLSB
		       << std::dec << std::setfill(' ')
		       << std::endl;
	    elIndx++;
	    if(elIndx%7==0) {
	      iel = 0;
	      ibx++;
	    }
	  }
	  /////////////////////////////////////////////////////////////////

	  ///////////////////// Fill to econtarray ///////////////////////
	  TPGBEDataformat::Trig24Data trdata[5][3]; //(3:trains+1:motherboard+1:external-trigger), 3:econts
	  for(int ilp=0;ilp<5;ilp++){
	    int nofecons = 0;
	    if(ilp==4)
	      nofecons = 1;
	    else if(ilp==3 or ilp==2)
	      nofecons = 2;
	    else
	      nofecons = 3;
	    for(int iecon=0;iecon<nofecons;iecon++){
	      moduleId = pck.packModId(zside, sector, ilp, det, iecon, selTC4, module);
	      if(ilp==4)
		trdata[ilp][iecon].setNofElinks( 7 );
	      else if(ilp==3)
		trdata[ilp][iecon].setNofElinks( ((iecon==0)?4:3) );
	      else if(ilp==2)
		trdata[ilp][iecon].setNofElinks( ((iecon==0)?4:2) );
	      else if(ilp==1)
		trdata[ilp][iecon].setNofElinks( 2 );
	      else
		trdata[ilp][iecon].setNofElinks( ((iecon==0)?3:2) );
	      trdata[ilp][iecon].setNofUnpkWords(8);
	      if((ilp==2 and iecon==0) or ilp==3) trdata[ilp][iecon].setNofUnpkWords(1,8);
	      for(uint32_t ib=0;ib<7;ib++){
		for(uint32_t iel=0;iel<trdata[ilp][iecon].getNofElinks();iel++){
		  if(iecon==0)
		    trdata[ilp][iecon].setElink(ib, iel, elpckt[ilp][ib][iel]) ;
		  else if(iecon==1)
		    if(ilp==1)
		      trdata[ilp][iecon].setElink(ib, iel, elpckt[ilp][ib][iel+trdata[ilp][0].getNofElinks()+1]) ;
		    else
		      trdata[ilp][iecon].setElink(ib, iel, elpckt[ilp][ib][iel+trdata[ilp][0].getNofElinks()]) ;
		  else
		    if(ilp==1)
		      trdata[ilp][iecon].setElink(ib, iel, elpckt[ilp][ib][iel+trdata[ilp][0].getNofElinks()+trdata[ilp][1].getNofElinks()+1]) ;
		    else
		      trdata[ilp][iecon].setElink(ib, iel, elpckt[ilp][ib][iel+trdata[ilp][0].getNofElinks()+trdata[ilp][1].getNofElinks()]) ;
		}
		for(uint8_t iw=0;iw<trdata[ilp][iecon].getNofUnpkWords();iw++) trdata[ilp][iecon].setUnpkWord(ib, iw, unpackedWord[ilp][iecon][ib][0][iw]) ;
		if((ilp==2 and iecon==0) or ilp==3) for(uint8_t iw=0;iw<trdata[ilp][iecon].getNofUnpkWords(1);iw++) trdata[ilp][iecon].setUnpkWord(ib, 1, iw, unpackedWord[ilp][iecon][ib][1][iw]) ;
	      }//nof bxs
	      trdata[ilp][iecon].setSlinkBx(eoe->bxId());
	      econtarray[eventId].push_back( std::make_pair(moduleId,trdata[ilp][iecon]) );
	    }//nof ECONTs
	  }//nof lpGBTs
	  ///////////////////////////////////////////////////////////////
	  
	  if(maxCAFESeps>=14){
	    //////////// TC proc block 10-13  /////////////////////
	    uint32_t TCProcWord[4][4][7][8]; //4:TDAQs, 4:columns, 7:bxs, 8:words
	    for(int iblk=12;iblk<=13;iblk++){ //we start with 12/13 as it corresponds to LD1/2 trains
	      unpkIndx = 0;
	      iunpkw = 0;
	      ibx = 0;
	      int itdaq = (iblk==12)?0:1 ; //0 and 1 for TDAQ 12/13
	      for(int iw = loc[iblk]+1; iw <= (loc[iblk]+size[iblk]) ; iw++ ){
		uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
		uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
		uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
		uint32_t col3 = p64[iw] & 0xFFFF ;
		TCProcWord[itdaq][0][ibx][iunpkw] = col3; //0-15 bits
		TCProcWord[itdaq][1][ibx][iunpkw] = col2; //16-31 bits
		TCProcWord[itdaq][2][ibx][iunpkw] = col1; //32-47 bits
		TCProcWord[itdaq][3][ibx][iunpkw] = col0; //48-63 bits
		iunpkw++;
		if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
		  std::cout<<"iblock: " << iblk <<", iloc: "<< iw 
			   << std::hex
			   <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
			   <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
			   <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
			   <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
			   << std::dec << std::setfill(' ')
			   << std::endl;		
		unpkIndx++;
		if(unpkIndx%8==0) {ibx++; iunpkw=0;}
	      }//loop over words for 7 bxs
	    }//loop over TC proc TDAQ blocks
	  
	    for(int iblk=10;iblk<=11;iblk++){ //Next we read 10/11 as it corresponds to LD3/MB1 trains
	      unpkIndx = 0;
	      iunpkw = 0;
	      ibx = 0;
	      int itdaq = (iblk==10)?2:3 ; //2 and 3 for TDAQ 10/11
	      for(int iw = loc[iblk]+1; iw <= (loc[iblk]+size[iblk]) ; iw++ ){
		uint32_t col0 = (p64[iw] >> (32+16)) & 0xFFFF ;
		uint32_t col1 = (p64[iw] >> 32) & 0xFFFF ;
		uint32_t col2 = (p64[iw] >> (32-16)) & 0xFFFF ;
		uint32_t col3 = p64[iw] & 0xFFFF ;
		TCProcWord[itdaq][0][ibx][iunpkw] = col3; //0-15 bits
		TCProcWord[itdaq][1][ibx][iunpkw] = col2; //16-31 bits
		TCProcWord[itdaq][2][ibx][iunpkw] = col1; //32-47 bits
		TCProcWord[itdaq][3][ibx][iunpkw] = col0; //48-63 bits
		iunpkw++;
		if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
		  std::cout<<"iblock: " << iblk <<", iloc: "<< iw 
			   << std::hex
			   <<", col0 : 0x" << std::setfill('0') << std::setw(4) << col0 <<", "
			   <<", col1 : 0x" << std::setfill('0') << std::setw(4) << col1 <<", "
			   <<", col2 : 0x" << std::setfill('0') << std::setw(4) << col2 <<", "
			   <<", col3 : 0x" << std::setfill('0') << std::setw(4) << col3 <<", "
			   << std::dec << std::setfill(' ')
			   << std::endl;	
		unpkIndx++;
		if(unpkIndx%8==0) {ibx++; iunpkw=0;}
	      }//loop over words for 7 bxs
	    }//loop over TC proc TDAQ blocks
	  
	    ///////////////////// Fill tcprocarray ///////////////////////
	    TPGBEDataformat::TrigTCProcData tcprocdata[4][3]; //4:lpGBT, 3:econT
	    for(int ilp=0;ilp<4;ilp++){
	      unsigned maxecons = ilp<2 ? 3:2;
	      for(int iecon=0;iecon< maxecons;iecon++){//2 ECON-Ts in trains 3 and 4
		moduleId = pck.packModId(zside, sector, ilp, det, iecon, selTC4, module);
		int theTDAQEntry = tcprocdata[ilp][iecon].getTDAQEntry(ilp,iecon);
		std::map<int, std::vector<std::pair<int, int>>> theWordAndColPerBin = tcprocdata[ilp][iecon].getWordAndColPerBin(ilp,iecon);
		for(uint32_t ib=0;ib<7;ib++){
		  for(uint32_t ibin=0;ibin<theWordAndColPerBin.size();ibin++){ 
		    uint32_t iinst=0;
		    for(const auto& thePair : theWordAndColPerBin[ibin]){
		      tcprocdata[ilp][iecon].setUnpkWord(ib, ibin, iinst, TCProcWord[theTDAQEntry][thePair.second][ib][thePair.first]) ;
		      iinst++;
		    }
		  }
		}//nof bxs
		tcprocarray[eventId].push_back( std::make_pair(moduleId,tcprocdata[ilp][iecon]) );
	      }//nof ECONTs
	    }//nof lpGBTs
	    ///////////////////////////////////////////////////////////////
	  }
	  
	}//MinMaxEvent
	if((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)) std::cout<<"========= End of event : "<< nEvents << "============="<<  std::endl;
	ievent++;
      }//loop event

    }//while read r
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONTReader::terminate(){
    _fileReader.close();
    if(!r) delete r;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}

#endif
