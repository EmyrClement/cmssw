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
    void getEvents(uint64_t&, uint64_t&, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>&, std::vector<uint64_t>&);
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
    int find_ffsep_word(const Hgcal10gLinkReceiver::RecordRunning *rEvent, int& isMSB, int n, int ffsep_word_loc = -1) {
      const uint64_t *p64(((const uint64_t*)rEvent)+1);      
      if (ffsep_word_loc > 0) {
	if (is_ffsep_word(p64[ffsep_word_loc], isMSB)) {
	  return ffsep_word_loc;
	}else
	  return ffsep_word_loc;
      } 
      else {
	if(scanMode and eventId==inspectEvent) std::cout << "hgc_roc::ECONDReader::find_ffsep_word(): missed" << std::endl;
	;
      }
 
      int cafe_counter = 0;
      int cafe_word_idx = -1;
      for(uint32_t i(0);i<rEvent->payloadLength();i++){
	const uint64_t word = p64[i];
	if (is_ffsep_word(word,isMSB)) { // if word == 0xfeca
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
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::init(uint32_t relayNumber, uint32_t runNumber, uint32_t linkNumber){
    //Make the buffer space for the records
    r = new Hgcal10gLinkReceiver::RecordT<4095>;
    _fileReader.setDirectory(std::string("dat/Relay")+std::to_string(relayNumber));
    _fileReader.openRun(runNumber,linkNumber);  
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONDReader::getEvents(uint64_t& minEventDAQ, uint64_t& maxEventDAQ, std::map<uint64_t,std::vector<std::pair<uint32_t,TPGFEDataformat::HalfHgcrocData>>>& hrocarray, std::vector<uint64_t>& events){
    
    //Set up specific records to interpet the formats
    const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
    const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    uint64_t nEvents = 0;
    
    std::map<uint32_t,TPGFEConfiguration::ConfigEconD>& econDPar = configs.getEconDPara();
    uint32_t idx = pck.packModId(zside, sector, link, det, econt, selTC4, module);
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    std::string modName = modNameMap.at(std::make_tuple(det, selTC4, module));
    const std::map<std::pair<std::string,std::tuple<uint32_t,uint32_t,uint32_t>>,uint32_t>&seqToRocpin = (pck.getDetType()==0)?configs.getSiSeqToROCpin():configs.getSciSeqToROCpin();
    
    int nRx = int(econDPar[idx].getNeRx());
    
    //Use the fileReader to read the records
    while(_fileReader.read(r)) {
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
	eventId = boe->eventId();
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){ 
	  event_dump(rEvent);
	  rEvent->RecordHeader::print();
	  boe->print();
	  eoe->print();
	  std::cout<<"Payload length : "<< rEvent->payloadLength() << std::endl;
	}
	
	int isMSB[6];
	int ffsep_word_loc[6];
	bool is_ff_sep_notfound = false;
	for(int iloc=0;iloc<6;iloc++){
	  isMSB[iloc] = -1;
	  ffsep_word_loc[iloc] = find_ffsep_word(rEvent, isMSB[iloc], iloc+1);
	  if(isMSB[iloc] == -1) is_ff_sep_notfound = true;
	  if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){ 
	    if(iloc==0)
	      std::cout << "\nEvent : " << nEvents << ": " <<iloc+1 <<"-st 0xffff 0xffff separator word at location : " << ffsep_word_loc[iloc] << " word and isMSB : " << isMSB[iloc] << std::endl;
	    else if(iloc==1)
	      std::cout << "Event : " << nEvents << ": " <<iloc+1 <<"-nd 0xffff 0xffff separator word at location : " << ffsep_word_loc[iloc] << " word and isMSB : " << isMSB[iloc] << std::endl;
	    else if(iloc==2)
	      std::cout << "Event : " << nEvents << ": " <<iloc+1 <<"-rd 0xffff 0xffff separator word at location : " << ffsep_word_loc[iloc] << " word and isMSB : " << isMSB[iloc] << std::endl;
	    else
	      std::cout << "Event : " << nEvents << ": " <<iloc+1 <<"-th 0xffff 0xffff separator word at location : " << ffsep_word_loc[iloc] << " word and isMSB : " << isMSB[iloc] << std::endl;
	  }
	}
	
	int expctd_empty_word_loc = ffsep_word_loc[5]+19; //19 is valid only for full ROCs, in case of partials this will be different, but could be calculated from mapping (here floor(37/2)+1 = 19)
	int empty_isMSB = -1;
	bool hasfound_emptry_word = found_empty_word(rEvent, empty_isMSB, expctd_empty_word_loc) ;
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent))
	  std::cout << "Event : " << nEvents << ": 0x0000 0x0000 closing separator word at location : " << expctd_empty_word_loc << " hasfound_emptry_word : "<< hasfound_emptry_word << ",  and isMSB : " << empty_isMSB << std::endl;
      
	if(!hasfound_emptry_word) continue;
	if(is_ff_sep_notfound) continue;
	//if(nEvents>maxEvent) continue;
      
	int ichip = 0;
	const uint64_t *p64(((const uint64_t*)rEvent)+1);
      
	uint32_t wordH = ((p64[ffsep_word_loc[0]-1] >> 32) & 0xFFFFFFFF) ;
	int daq_header_payload_length = int(((wordH>>14) & 0x1FF));
	daq_header_payload_length -= 1 ; //1 for DAQ header
	if ((nEvents < nShowEvents) or (scanMode and boe->eventId()==inspectEvent)){
	  std::cout << std::hex << std::setfill('0');
	  std::cout << "Event : " << nEvents << ", WordH : 0x" << wordH  ;
	  std::cout << std::dec << std::setfill(' ');
	  std::cout<<", Size in ECON-D header : " << daq_header_payload_length << std::endl;
	}
      
	int record_header_sizeinfo = int(rEvent->payloadLength()); //in 64-bit format
	int reduced_record_header_sizeinfo = record_header_sizeinfo - (4+2+1) ; //exclude 4 64-bit words for slink header/trailer + 2 DAQ header + 1 for DAQ trailer
	int record_header_sizeinfo_32bit = 2*reduced_record_header_sizeinfo ; //change to 32-bit format

	if(record_header_sizeinfo_32bit!=daq_header_payload_length) continue;
      
	if(boe->eventId()>=minEventDAQ and boe->eventId()<maxEventDAQ){

	
	  for(int iloc=0;iloc<nRx;iloc++){ //For full LD/HD modules : 6/12 respectively
	    uint32_t max_sep_word = (iloc==5) ? expctd_empty_word_loc : ffsep_word_loc[iloc+1] ;
	    bool ishalf = (iloc%2==0) ? false : true ;
	    if(isMSB[iloc]==1) max_sep_word++;
	    int index = 0;
	    uint32_t ch = 0;
	    uint32_t rocpin = 0;
	    int adcmL = 0, adcL = 0, toaL = 0, totL = 0;
	    int adcmM = 0, adcM = 0, toaM = 0, totM = 0;
	    rocn = uint32_t(ichip);
	    half = uint32_t(ishalf);
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
		if(trigflagM>=0x2)
		  chdata[ch].setTot(uint16_t(totM),uint16_t(trigflagM));
		else if(trigflagM==0)
		  chdata[ch].setAdc(uint16_t(adcM),uint16_t(trigflagM));
		else if(trigflagM==0x1)
		  chdata[ch].setAdc(0,uint16_t(trigflagM)); //uint16_t(adcM)
		else
		  chdata[ch].setZero();
	      }
	      if(seqToRocpin.find(std::make_pair(modName,std::make_tuple(rocn,half,seqL)))!=seqToRocpin.end()){
		uint32_t rocpin = seqToRocpin.at(std::make_pair(modName,std::make_tuple(rocn,half,seqL))) ;
		ch = rocpin%36 ; //36 for halfroc
		if(trigflagL>=0x2)
		  chdata[ch].setTot(uint16_t(totL),uint16_t(trigflagL));
		else if(trigflagL==0)
		  chdata[ch].setAdc(uint16_t(adcL),uint16_t(trigflagL));
		else if(trigflagL==0x1)
		  chdata[ch].setAdc(0,uint16_t(trigflagL)); //uint16_t(adcL)
		else
		  chdata[ch].setZero();
	      }
	      index++;
	    }//end of loop for 37 seq channels
	    TPGFEDataformat::HalfHgcrocData hrocdata;
	    hrocdata.setChannels(chdata);
	    pck.setZero();
	    hrocarray[boe->eventId()].push_back(std::make_pair(pck.packRocId(zside, sector, link, det, econt, selTC4, module, rocn, half),hrocdata));
	    //std::cout<<std::endl;
	    if(iloc%2==1)ichip++;
	  }
	  events.push_back(boe->eventId());
	}
	//Increment event counter
	nEvents++;	
      }
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
    
    ECONTReader(TPGFEConfiguration::Configuration& cfgs) : configs(cfgs), r(0x0), scanMode(false), inspectEvent(0), nShowEvents(0), eventId(0), moduleId(0) { initId();}
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
    
    void init(uint32_t, uint32_t, uint32_t);
    void getEvents(uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<TPGFEDataformat::TcModulePacket>>& econtarray, std::vector<uint64_t>& events){
      moduleId = pck.packModId(zside, sector, link, det, econt, selTC4, module);
      const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
      if(outputType==TPGFEDataformat::BestC)
	getEventsBC(minEventTrig, maxEventTrig, econtarray, events);
      else if(outputType==TPGFEDataformat::STC4A or outputType==TPGFEDataformat::STC4B or outputType==TPGFEDataformat::STC16)
	getEventsSTC(minEventTrig, maxEventTrig, econtarray, events);
      else{
	std::vector<TPGFEDataformat::TcModulePacket> dummy1;
	econtarray[0] = dummy1;
    }

    }
    void getEventsBC(uint64_t&, uint64_t&, std::map<uint64_t,std::vector<TPGFEDataformat::TcModulePacket>>&, std::vector<uint64_t>&);
    void getEventsSTC(uint64_t&, uint64_t&, std::map<uint64_t,std::vector<TPGFEDataformat::TcModulePacket>>&, std::vector<uint64_t>&);
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
    
    bool is_cafe_word(const uint32_t word) {
      if(((word >> 8) & 0xFFFFFF) == 0xFECAFE)
	return true; //0xFECAFE                                                      
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
	;//std::cout << "oops" << std::endl;
      }
 
      int cafe_counter = 0;
      int cafe_word_idx = -1;
      for(unsigned i(0);i<rEvent->payloadLength();i++){
	const uint32_t word = p64[i];
	if (is_cafe_word(word)) { // if word == 0xfecafe
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
    uint32_t trig_linkNumber = TMath::FloorNint((linkNumber-1)/2);
    r = new Hgcal10gLinkReceiver::RecordT<4095>;
    _fileReader.setDirectory(std::string("dat/Relay")+std::to_string(relayNumber));
    _fileReader.openRun(runNumber,trig_linkNumber);
    isMSB = (linkNumber%2==0) ? true : false;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONTReader::getEventsSTC(uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<TPGFEDataformat::TcModulePacket>>& econtarray, std::vector<uint64_t>& events){

    //Set up specific records to interpet the formats
    const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
    const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    
    TPGFEDataformat::TcRawDataPacket edata;
    std::vector<TPGFEDataformat::TcModulePacket> edataarray;
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    const uint32_t maxSTcs = (pck.getDetType()==0)?configs.getSiModSTClist().at(modName).size():configs.getSciModSTClist().at(modName).size();
    
    uint64_t nEvents = 0;
    uint64_t prevEvent = 0;
    uint32_t prevSequence = 0;
    uint32_t packet[4];
    uint32_t packet_counter;
    uint32_t packet_locations[maxSTcs];
    uint64_t packet_energies[maxSTcs];
  
    uint16_t daq_data[5];                              
    uint16_t daq_nbx[5];                               
    uint16_t size_in_cafe[5];                          
    uint32_t bx_raw[2][15];                        
    int refbxindex = 7;
    
    //Use the fileReader to read the records
    while(_fileReader.read(r)) {
      //Check the state of the record and print the record accordingly
      if( r->state()==Hgcal10gLinkReceiver::FsmState::Starting){
	if(!(rStart->valid())){
	  std::cerr << " FsmState::Starting validadity fails : rStart->valid() " << rStart->valid() << std::endl;
	  // rStart->print();
	  // std::cout << std::endl;
	  continue;
	}
      }
    
      else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping){
	if(!(rStop->valid())){
	  std::cerr << " FsmState::Stopping validadity fails : rStop->valid() " << rStop->valid() << std::endl;
	  // rStop->print();
	  // std::cout << std::endl;
	  continue;
	}
      }
      //Else we have an event record 
      else{
      
	const Hgcal10gLinkReceiver::SlinkBoe *boe = rEvent->slinkBoe();      
	const Hgcal10gLinkReceiver::SlinkEoe *eoe = rEvent->slinkEoe();
	eventId = boe->eventId();
	edata.reset();
	edata.setTBM(outputType, 0, 0);
	edataarray.clear();
	// if (nEvents < 1){ 
	//  	event_dump(rEvent);
	// 	rEvent->RecordHeader::print();
	// 	boe->print();
	// 	eoe->print();
	// }

	if(boe->boeHeader()!=boe->BoePattern) continue;
	if(eoe->eoeHeader()!=eoe->EoePattern) continue;
	uint16_t l1atype = boe->l1aType();      
	//if(l1atype==0) continue;
      
	const uint64_t *p64(((const uint64_t*)rEvent)+1);

	// edata.event = nEvents;
	// edata.eventId = eventId;
	// edata.bxId = eoe->bxId();
	// edata.sequenceId = rEvent->RecordHeader::sequenceCounter(); 

	if((Abs64(eventId,prevEvent) != Abs32(rEvent->RecordHeader::sequenceCounter(), prevSequence)) and Abs64(eventId,prevEvent)>=2){
	  prevEvent = eventId;
	  prevSequence = rEvent->RecordHeader::sequenceCounter(); 
	  continue;
	}
	prevEvent = eventId;
	prevSequence = rEvent->RecordHeader::sequenceCounter(); 

	daq_data[0] = p64[2] & 0xF;
	daq_nbx[0] = p64[2]>>4 & 0x7;
	uint64_t daq0_event_size = (2*daq_nbx[0] + 1)*daq_data[0];
	int first_cafe_word_loc = find_cafe_word(rEvent, 1);
	size_in_cafe[0] = p64[first_cafe_word_loc] & 0xFF;      
      
	daq_data[1] = p64[2]>>7 & 0xF;
	daq_nbx[1] = p64[2]>>11 & 0x7;
	uint64_t daq1_event_size = (2*daq_nbx[1] + 1)*daq_data[1];
	int second_cafe_word_loc = find_cafe_word(rEvent, 2);
	size_in_cafe[1] = p64[second_cafe_word_loc] & 0xFF;
	
	int sixth_cafe_word_loc = find_cafe_word(rEvent, 6);
	
	if(sixth_cafe_word_loc!=0) continue;
	if(first_cafe_word_loc != 3) continue;
	if(daq0_event_size != size_in_cafe[0]) continue;	
	if(daq1_event_size != size_in_cafe[1]) continue;	
	
	if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
	  event_dump(rEvent);
	  rEvent->RecordHeader::print();
	  boe->print();
	  eoe->print();
	}
	//if(nEvents>maxEvent) continue;      
	if(eventId>=minEventTrig and eventId<maxEventTrig){
	  
	  int bx_index = -1.0*int(daq_nbx[0]);
	  const int maxnbx = (2*daq_nbx[0] + 1);
	  uint64_t energy_raw[2][maxnbx][maxSTcs];
	  for(int iect=0;iect<2;iect++)
	    for(int ibx=0;ibx<maxnbx;ibx++){
	      bx_raw[iect][ibx] = 0;
	      for(int istc=0;istc<maxSTcs;istc++){
		energy_raw[iect][ibx][istc] = 0;
		// edata.energy_raw[iect][istc] = 0;
		// edata.loc_raw[iect][istc] = 0;
	      }
	    }

	  uint32_t word, bx_counter ;

	  for(unsigned i(first_cafe_word_loc+1);i<daq0_event_size+first_cafe_word_loc+1;i=i+4){

	    if(!isMSB){
	      
	      word = (p64[i] & 0xFFFFFFFF) ;
	      bx_counter = ( word >> 28 ) & 0xF;
	
	      packet[0] = (p64[i] & 0xFFFFFFFF) ;
	      packet[1] = (p64[i+1] & 0xFFFFFFFF) ;
	      packet[2] = (p64[i+2] & 0xFFFFFFFF) ; 
	      packet[3] = (p64[i+3] & 0xFFFFFFFF) ;
	
	      set_packet_locations(packet_locations, packet);
	      set_packet_energies(packet_energies, packet);
	
	      uint64_t totEnergy = 0;
	      for(int istc=0;istc<maxSTcs;istc++) totEnergy += packet_energies[istc] ;
	
	      if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
		std::cout<<" EventId : "<<eventId
			 << ", bx_index : " << (bx_index+int(daq_nbx[0])) 
			 <<", bx(LSB) :"  << bx_counter
			 <<", edata.bxId : "<<eoe->bxId() <<", modulo8 : "<< (eoe->bxId()%8)
			 <<", totEnergy : "<<totEnergy
			 <<std::endl;
		std::cout<<"E(LSB) : \t";
		for(int istc=0;istc<maxSTcs;istc++){
		  std::cout<<packet_energies[istc]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"L(LSB) : \t";
		for(int istc=0;istc<maxSTcs;istc++){
		  std::cout<<packet_locations[istc]<<" ";
		}
		std::cout<<std::endl;
	      }
	
	      bx_raw[0][(bx_index+int(daq_nbx[0]))] = bx_counter;
	      for(int istc=0;istc<maxSTcs;istc++){
		energy_raw[0][(bx_index+int(daq_nbx[0]))][istc] = packet_energies[istc];
	      }
	    }else{
	
	
	      word = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      bx_counter = ( word >> 28 ) & 0xF;
	
	      packet[0] = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      packet[1] = ((p64[i+1] >> 32) & 0xFFFFFFFF) ;
	      packet[2] = ((p64[i+2] >> 32) & 0xFFFFFFFF) ;
	      packet[3] = ((p64[i+3] >> 32) & 0xFFFFFFFF) ;
	
	      set_packet_locations(packet_locations, packet);
	      set_packet_energies(packet_energies, packet);

	      uint64_t totEnergy1 = 0;
	      for(int istc=0;istc<maxSTcs;istc++) totEnergy1 += packet_energies[istc] ;

	      if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
		std::cout<<" EventId : "<<eventId
			 << ", bx_index : " << (bx_index+int(daq_nbx[0])) 
			 <<", bx(MSB) :"  << bx_counter
			 <<", edata.bxId : "<<eoe->bxId() <<", modulo8 : "<< (eoe->bxId()%8)
			 <<", totEnergy : "<<totEnergy1
			 <<std::endl;
		std::cout<<"E(MSB) : \t";
		for(int istc=0;istc<maxSTcs;istc++){
		  std::cout<<packet_energies[istc]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"L(MSB) : \t";
		for(int istc=0;istc<maxSTcs;istc++){
		  std::cout<<packet_locations[istc]<<" ";
		}
		std::cout<<std::endl;
	      }
	
	      bx_raw[1][(bx_index+int(daq_nbx[0]))] = bx_counter;
	      for(int istc=0;istc<maxSTcs;istc++){
		energy_raw[1][(bx_index+int(daq_nbx[0]))][istc] = packet_energies[istc];
	      }
	    }//isMSB
	    
	    if((bx_index+int(daq_nbx[0]))==refbxindex){
	      for(int istc=0;istc<maxSTcs;istc++){
		//edata.push_back(TPGFEDataformat::TcRawData(outputType, , packet_energies[istc]));
		edata.setTcData(outputType, packet_locations[istc], packet_energies[istc]);
	      }
	    }
	    
	    bx_index++;
	  }//loop over unpacked

	  TPGFEDataformat::TcModulePacket modpacket = std::make_pair(moduleId,edata);
	  edataarray.push_back(modpacket);
    	  econtarray[eventId] = edataarray;
	  
	  //econtarray[eventId] = edata;
	  //if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) std::cout<<std::endl;
	}
	//Increment event counter
	nEvents++;

      }
    }//file reader

  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //void ECONTReader::getEventsBC(uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<std::pair<uint32_t,std::vector<TPGFEDataformat::TcRawData>>>>& econtarray, std::vector<uint64_t>& events){
  void ECONTReader::getEventsBC(uint64_t& minEventTrig, uint64_t& maxEventTrig, std::map<uint64_t,std::vector<TPGFEDataformat::TcModulePacket>>& econtarray, std::vector<uint64_t>& events){
    
    //Set up specific records to interpet the formats
    const Hgcal10gLinkReceiver::RecordStarting *rStart((Hgcal10gLinkReceiver::RecordStarting*)r);
    const Hgcal10gLinkReceiver::RecordStopping *rStop ((Hgcal10gLinkReceiver::RecordStopping*)r);
    const Hgcal10gLinkReceiver::RecordRunning  *rEvent((Hgcal10gLinkReceiver::RecordRunning*) r);
    
    TPGFEDataformat::TcRawDataPacket edata;
    std::vector<TPGFEDataformat::TcModulePacket> edataarray;
    
    const TPGFEDataformat::Type& outputType = configs.getEconTPara().at(moduleId).getOutType();
    const uint32_t maxBCTcs = configs.getEconTPara().at(moduleId).getBCType() ;
    
    const std::map<std::tuple<uint32_t,uint32_t,uint32_t>,std::string>& modNameMap = configs.getModIdxToName();
    const std::string& modName = modNameMap.at(std::make_tuple(pck.getDetType(),pck.getSelTC4(),pck.getModule()));
    const uint32_t maxNofTcs = (pck.getDetType()==0)?configs.getSiModTClist().at(modName).size():configs.getSciModTClist().at(modName).size();

    uint64_t nEvents = 0;
    uint64_t prevEvent = 0;
    uint32_t prevSequence = 0;
    uint32_t packet[4];
    uint32_t packet_counter;
    bool packet_tcs[maxNofTcs];
    uint64_t packet_energies[maxBCTcs];

    uint16_t daq_data[5];                              
    uint16_t daq_nbx[5];                               
    uint16_t size_in_cafe[5];                          
    uint32_t bx_raw[2][15];                        
    int refbxindex = 10;
  
    //Use the fileReader to read the records
    while(_fileReader.read(r)) {
      //Check the state of the record and print the record accordingly
      if( r->state()==Hgcal10gLinkReceiver::FsmState::Starting){
	if(!(rStart->valid())){
	  std::cerr << " FsmState::Starting validadity fails : rStart->valid() " << rStart->valid() << std::endl;
	  // rStart->print();
	  // std::cout << std::endl;
	  continue;
	}
      }
      
      else if(r->state()==Hgcal10gLinkReceiver::FsmState::Stopping){
	if(!(rStop->valid())){
	  std::cerr << " FsmState::Stopping validadity fails : rStop->valid() " << rStop->valid() << std::endl;
	  // rStop->print();
	  // std::cout << std::endl;
	  continue;
	}
      }
      //Else we have an event record 
      else{

	const Hgcal10gLinkReceiver::SlinkBoe *boe = rEvent->slinkBoe();      
	const Hgcal10gLinkReceiver::SlinkEoe *eoe = rEvent->slinkEoe();
	eventId = boe->eventId();
	//clean the array;
	edata.reset();
	edataarray.clear();

	// if (nEvents < 1){ 
	//  	event_dump(rEvent);
	// 	rEvent->RecordHeader::print();
	// 	boe->print();
	// 	eoe->print();
	// }
	
	
	if(boe->boeHeader()!=boe->BoePattern) continue;
	if(eoe->eoeHeader()!=eoe->EoePattern) continue;
	uint16_t l1atype = boe->l1aType();      
	//if(l1atype==0) continue;
	
	const uint64_t *p64(((const uint64_t*)rEvent)+1);
	
	if((Abs64(eventId,prevEvent) != Abs32(rEvent->RecordHeader::sequenceCounter(), prevSequence)) and Abs64(eventId,prevEvent)>=2){
	  prevEvent = eventId;
	  prevSequence = rEvent->RecordHeader::sequenceCounter(); 
	  continue;
	}
	prevEvent = eventId;
	prevSequence = rEvent->RecordHeader::sequenceCounter(); 
	
	daq_data[0] = p64[2] & 0xF;
	daq_nbx[0] = p64[2]>>4 & 0x7;
	uint64_t daq0_event_size = (2*daq_nbx[0] + 1)*daq_data[0];
	int first_cafe_word_loc = find_cafe_word(rEvent, 1);
	size_in_cafe[0] = p64[first_cafe_word_loc] & 0xFF;      
	
	daq_data[1] = p64[2]>>7 & 0xF;
	daq_nbx[1] = p64[2]>>11 & 0x7;
	uint64_t daq1_event_size = (2*daq_nbx[1] + 1)*daq_data[1];
	int second_cafe_word_loc = find_cafe_word(rEvent, 2);
	size_in_cafe[1] = p64[second_cafe_word_loc] & 0xFF;
	
	daq_data[2] = p64[2]>>14 & 0xF;
	daq_nbx[2] = p64[2]>>18 & 0x7;
	uint64_t daq2_event_size = (2*daq_nbx[2] + 1)*daq_data[2];
	int third_cafe_word_loc = find_cafe_word(rEvent, 3);
	size_in_cafe[2] = p64[third_cafe_word_loc] & 0xFF;
	
	int sixth_cafe_word_loc = find_cafe_word(rEvent, 6);
      
	if(sixth_cafe_word_loc!=0) continue;
	if(first_cafe_word_loc != 3) continue;
	if(daq0_event_size != size_in_cafe[0]) continue;	
	if(daq1_event_size != size_in_cafe[1]) continue;
	if(daq2_event_size != size_in_cafe[2]) continue;	
	
	if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
	  event_dump(rEvent);
	  rEvent->RecordHeader::print();
	  boe->print();
	  eoe->print();
	}
	
	if(eventId>=minEventTrig and eventId<maxEventTrig){
	  
    	  int bx_index = -1.0*int(daq_nbx[0]);
    	  const int maxnbx = (2*daq_nbx[0] + 1);
    	  uint64_t energy_raw[2][maxnbx][maxBCTcs];
    	  for(int iect=0;iect<2;iect++)
    	    for(int ibx=0;ibx<maxnbx;ibx++){
    	      bx_raw[iect][ibx] = 0;
    	      for(int istc=0;istc<maxBCTcs;istc++){
    	       	energy_raw[iect][ibx][istc] = 0;
    	      }
    	    }
	  
    	  for(unsigned i(first_cafe_word_loc+1);i<daq0_event_size+first_cafe_word_loc+1;i=i+4){

	    uint32_t word, bx_counter, modsum;
	    
	    if(!isMSB){
	      
	      word = (p64[i] & 0xFFFFFFFF) ;
	      bx_counter = ( word >> 28 ) & 0xF;
	      modsum = ( word >> 20 ) & 0xFF;
	    
	      packet[0] = (p64[i] & 0xFFFFFFFF) ;
	      packet[1] = (p64[i+1] & 0xFFFFFFFF) ;
	      packet[2] = (p64[i+2] & 0xFFFFFFFF) ; 
	      packet[3] = (p64[i+3] & 0xFFFFFFFF) ;
	      
	      set_packet_tcs_bc(packet_tcs, packet, maxNofTcs);
	      set_packet_energies_bc9(packet_energies, packet, maxBCTcs);
	      
	      if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
		std::cout<<" EventId : "<<eventId
			 <<", index : "<<(bx_index+int(daq_nbx[0]))
			 <<", bx(LSB) :"  << bx_counter 
			 <<", edata.bxId : "<<eoe->bxId() <<", modulo8 : "<< (eoe->bxId()%8)
			 <<", modsum : "<<modsum
    		       <<std::endl;
		std::cout<<"E(LSB) : \t";
		for(int itc=0;itc<maxBCTcs;itc++){
		  std::cout<<packet_energies[itc]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"L(LSB) : \t";
		for(int itc=0;itc<maxNofTcs;itc++){
		  if(packet_tcs[itc])
		    std::cout<<itc<<" ";
		}
		std::cout<<std::endl;
	      }
	      
	      bx_raw[0][(bx_index+int(daq_nbx[0]))] = bx_counter;
	      for(int itc=0;itc<maxBCTcs;itc++){	  
		energy_raw[0][(bx_index+int(daq_nbx[0]))][itc] = packet_energies[itc];
	      }
	      
	    }else{
	      
	      word = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      bx_counter = ( word >> 28 ) & 0xF;
	      modsum = ( word >> 20 ) & 0xFF;
	      
	      packet[0] = ((p64[i] >> 32) & 0xFFFFFFFF) ;
	      packet[1] = ((p64[i+1] >> 32) & 0xFFFFFFFF) ;
	      packet[2] = ((p64[i+2] >> 32) & 0xFFFFFFFF) ;
	      packet[3] = ((p64[i+3] >> 32) & 0xFFFFFFFF) ;
	      
	      set_packet_tcs_bc(packet_tcs, packet, maxNofTcs);
	      set_packet_energies_bc9(packet_energies, packet, maxBCTcs);
	      
	      if ((nEvents < nShowEvents) or (scanMode and eventId==inspectEvent)) {
		std::cout<<" EventId : "<<eventId
			 <<", index : "<<(bx_index+int(daq_nbx[0]))
			 <<", bx(MSB) :"  << bx_counter
			 <<", edata.bxId : "<<eoe->bxId() <<", modulo8 : "<< (eoe->bxId()%8)
			 <<", modsum : "<<modsum
			 <<std::endl;
		std::cout<<"E(MSB) : \t";
		for(int itc=0;itc<maxBCTcs;itc++){
		  std::cout<<packet_energies[itc]<<" ";
		}
		std::cout<<std::endl;
		std::cout<<"L(MSB) : \t";
		for(int itc=0;itc<maxNofTcs;itc++){
		  if(packet_tcs[itc])
		    std::cout<<itc<<" ";
		}
		std::cout<<std::endl;
	      }

	      bx_raw[1][(bx_index+int(daq_nbx[0]))] = bx_counter;
	      for(int itc=0;itc<maxBCTcs;itc++){	  
		energy_raw[1][(bx_index+int(daq_nbx[0]))][itc] = packet_energies[itc];
	      }
	    }
	    int itrigcell = 0;
	    if((bx_index+int(daq_nbx[0]))==refbxindex){
	      //edata.push_back(TPGFEDataformat::TcRawData(modsum));
	      edata.setTBM(outputType, 0, modsum) ;
	      for(int itc=0;itc<maxNofTcs;itc++)
    		if(packet_tcs[itc]){
		  //edata.push_back(TPGFEDataformat::TcRawData(outputType, itc, packet_energies[itrigcell++]));
		  //std::cout << "itc: " << itc << ", maxNofTcs: "<<maxNofTcs<<", e: " << packet_energies[itrigcell] << std::endl;
		  edata.setTcData(outputType, itc, packet_energies[itrigcell++]);
		}
	    }
	    
    	    //Itrigcell = 0;
    	    // if((bx_index+int(daq_nbx[0]))==refbxindex){
    	    //   edata.modsum[1] = modsum1;
    	    //   for(int itc=0;itc<9;itc++)
    	    // 	edata.energy_raw[1][itc] = packet_energies[itc];
    	    //   for(int itc=0;itc<48;itc++)
    	    // 	if(packet_tcs[itc]) edata.loc_raw[1][itrigcell++] = itc;
    	    // }
	
    	    bx_index++;
    	  }
	  TPGFEDataformat::TcModulePacket modpacket = std::make_pair(moduleId,edata);
	  edataarray.push_back(modpacket);
    	  econtarray[eventId] = edataarray;
	}
	//Increment event counter
	nEvents++;
      }
    }//file reader
    
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ECONTReader::terminate(){
    _fileReader.close();
    if(!r) delete r;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
}

#endif
