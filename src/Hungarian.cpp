
#include "Hungarian.h"


void
Hungarian::HungarianAlgo(BipartiteGraph& _bg, vector<VID>& _S, vector<VID>& _T,
                        vector<VID>& _N, vector<EID>& _EG, vector<EID>& _M){

  int MaxLoop1 = _bg.GetNumTasks()*2;
  int MaxLoop2 = _bg.GetNumTasks()*_bg.GetNumTasks()*_bg.GetNumTasks();
  int cnt1 = 0;
  while(cnt1++ < MaxLoop1){
    //DisplayData(1,1,1,1,1);
    if(IsPerfect(_bg)) break;
    
    VID u = PickFreeAgent(_bg);
    InitNewRoot(_bg, u, _S, _T);
    
    int cnt2 = 0;
    while(cnt2++ < MaxLoop2){
      RefreshBG(_bg, _S, _T, _N, _EG, _M);
      
      //if need relabel, update labels
      if(NeedReLabel(_T, _N)){
	UpdateLabels(_bg);
      }
      RefreshBG(_bg, _S, _T, _N, _EG, _M);

      //if not need relabel
      if(!NeedReLabel(_T, _N)){
	VID y = PickAvailableTask(_T, _N);
 	
	//if NOT need relabel and task NOT matched, then found augmenting path
	if(!bg.GetTask(y)->GetMatched()){
	  vector<EID> _path = BFSAugmentingPath(_bg, u, y);
	  AugmentPath(_bg, _path);
	  RefreshBG(_bg, _S, _T, _N, _EG, _M);

	  break; //break innter while
        }
  	//if not need relabel and task is *matched*, extend tree
        else{
	  ExtendTree(_bg, y, _S, _T);
        }
      }//if(!NeedReLabel)
    } //while(cnt2)
    if(cnt2 >= MaxLoop2){
      cerr<<"I suspect the correctness of Hungarian inner loop. Given up."<<endl;
      exit(-1);
    }
  }//while(cnt1)
  if(cnt1 >= MaxLoop1){
    cout<<"I suspect the correctness of Hungarian external loop. Given up."<<endl;
    exit(-1);
  }
  
  cout<<endl<<"*****************************************************"<<endl;
  cout<<endl<<"Your queried assignment problem is:"<<endl;
  DisplayMatrix(*(_bg.GetMatrix()));
  cout<<"The matching configuration is:"<<endl;
  DisplayConfig(_bg);

  cout<<"@_@ Kuhn-Munkres Hungarian Perfect Matching:"<<endl;
  DisplayData(_M);
  cout<<"Optimization result (weights sum): "<<this->OptimalValue(_bg, _M);
  cout<<endl;
  cout<<endl<<"*****************************************************"<<endl;
  cout<<endl;

}



void
Hungarian::HungarianAlgo(vector<EID>& matchtable){

  this->HungarianAlgo(this->bg, this->S, this->T, this->N,
		      this->EG, this->M);
	matchtable = M;

}



void
Hungarian::InitNewRoot(BipartiteGraph& _bg, VID root, vector<VID>& _S, vector<VID>& _T){
 
  //init the sets of S and T
  _S.clear();
  _T.clear();
  _S.push_back(root);
  
  //clear history
  for(size_t i=0; i<_bg.GetNumAgents(); i++)
    _bg.GetAgent(i)->SetColored(false);
  for(size_t j=0; j<_bg.GetNumTasks(); j++)
    _bg.GetTask(j)->SetColored(false);

  //color the root
  _bg.GetAgent(root)->SetColored(true);

}


void
Hungarian::InitNewRoot(VID root){

  this->InitNewRoot(this->bg, root, this->S, this->T);

}


void 
Hungarian::RefreshBG(BipartiteGraph& _bg, 
		const vector<VID>& _S, const vector<VID>& _T, 
		vector<VID>& _N, vector<EID>& _EG, vector<EID>& _M){
   
  //check feasibility, if not, warning
  if(!_bg.CheckFeasibility())
    cout<<"Warning: Bipartite Graph is not feasible!"<<endl;

  //check admissible edges, update admissible flags for edges, and update set EG
  _EG.clear();
  for(unsigned int i=0; i<_bg.GetNumAgents(); i++)
    for(unsigned int j=0; j<_bg.GetNumTasks(); j++){
      if( fabs(_bg.GetMatrix(i,j)->GetWeight() - (_bg.GetAgent(i)->GetLabel() + _bg.GetTask(j)->GetLabel())) <= DOUBLE_EPSILON ){
	_bg.GetMatrix(i,j)->SetAdmissibleFlag(true);
	_EG.push_back(EID(i,j));
      }
      else
	_bg.GetMatrix(i,j)->SetAdmissibleFlag(false);
    } //for j

  //check matched edges, update matched flags for vertices, and update set M
  _M.clear();
  size_t cnt = 0;
  for(unsigned int i=0; i<_bg.GetNumAgents(); i++)
    _bg.GetAgent(i)->SetMatched(false);
  for(unsigned int j=0; j<_bg.GetNumTasks(); j++)
    _bg.GetTask(j)->SetMatched(false);
  for(unsigned int i=0; i<_bg.GetNumAgents(); i++)
    for(unsigned int j=0; j<_bg.GetNumTasks(); j++)
      if(_bg.GetMatrix(i,j)->GetMatchedFlag()){
	_bg.GetAgent(i)->SetMatched(true);
	_bg.GetTask(j)->SetMatched(true); 
	_M.push_back(EID(i, j));
 	cnt++;
      }//if

  //update the variable in BG
  _bg.SetNumMatched(cnt); 

  //update set N
  _N.clear();
  _N = this->FindNeighbors(_EG, _S); 
 
}

void 
Hungarian::RefreshBG(void){

  this->RefreshBG(this->bg, this->S, this->T, this->N, this->EG, this->M);

}



void 
Hungarian::ExtendTree(BipartiteGraph& _bg, VID x, vector<VID>& _S, vector<VID>& _T){
  //if task x matched some agent
  int vid_agent = -1;
  if(_bg.GetTask(x)->GetMatched()){
    for(unsigned int i=0; i<_bg.GetNumAgents(); i++)
      if(_bg.GetMatrix(i, x)->GetMatchedFlag()){
	vid_agent = i;
	break;
      }

    if(-1 == vid_agent){
      cerr<<"Error: Matched task vertex could not find its matching agent. Stopped."<<endl;
      exit(-1);
    }

    _S.push_back(vid_agent);
    _bg.GetAgent(vid_agent)->SetColored(true);
    _T.push_back(x);
    _bg.GetTask(x)->SetColored(true);
    //EG.push_back(EID(vid_agent,x));
  }//if

}


void 
Hungarian::ExtendTree(VID x){

  this->ExtendTree(this->bg, x, this->S, this->T);

}


//AugmentPath should be followed by RefreshBG(), since after this, config changes.
void 
Hungarian::AugmentPath(BipartiteGraph& _bg, vector<EID>& _path){
  for(vector<EID>::iterator itr = _path.begin(); itr != _path.end(); itr++){
    if(_bg.GetMatrix((*itr).first, (*itr).second)->GetMatchedFlag())
      _bg.GetMatrix((*itr).first, (*itr).second)->SetMatchedFlag(false);
    else
      _bg.GetMatrix((*itr).first, (*itr).second)->SetMatchedFlag(true);
  }

}

 
void 
Hungarian::AugmentPath(vector<EID>& _path){

  this->AugmentPath(this->bg, _path);

}




vector<EID> 
Hungarian::BFSAugmentingPath(BipartiteGraph& _bg, VID x, VID y){
  size_t tasks_size = _bg.GetNumTasks();
  size_t agents_size = _bg.GetNumAgents();
  bool found = false;
  vector<EID> aug_path;
  deque<Vertex> dq;

  //before alternating tree, clear all paths
  for(unsigned int i=0; i<tasks_size; i++) 
    _bg.GetTask(i)->path.clear();
  for(unsigned int j=0; j<agents_size; j++)
    _bg.GetAgent(j)->path.clear();

  dq.push_back(*_bg.GetAgent(x));
  //cout<<"dq.front: "<<dq.front().GetVID()<<endl;

  //loop for searching a path, using queue.
  while(dq.size() && !found){
    Vertex v = dq.front();
    dq.pop_front();
    if(v.GetObj() == "TASK"){
      size_t agent_index = 0;
      for(agent_index=0; agent_index<agents_size; agent_index++)
        if(_bg.GetMatrix(agent_index, v.GetVID())->GetMatchedFlag() && !_bg.GetAgent(agent_index)->path.size()){
	  _bg.GetAgent(agent_index)->path = v.path;
	  _bg.GetAgent(agent_index)->path.push_back(EID(agent_index, v.GetVID()));
 	  dq.push_back(*_bg.GetAgent(agent_index));
 	  break;          //only ONE possible matched edge for ONE vertex.
        }

      if(agent_index ==  agents_size){
	/*if(_Verbose)
	  cout<<"Warning: I found no matched edge during back-tracking from task vertex: "<< v.GetVID()<<endl
	      <<"         It's possible there exists new other satisfied task vertices. Continue."<<endl;*/
        //exit(-1);
      }
    }//if "TASK"

    else if(v.GetObj() == "AGENT"){
      size_t task_index = 0;
      for(task_index=0; task_index<tasks_size; task_index++)
	if(!_bg.GetMatrix(v.GetVID(), task_index)->GetMatchedFlag() 
		&& _bg.GetMatrix(v.GetVID(), task_index)->GetAdmissibleFlag() 
		&& !_bg.GetTask(task_index)->path.size()){
	  _bg.GetTask(task_index)->path = v.path;
	  _bg.GetTask(task_index)->path.push_back(EID(v.GetVID(), task_index));
	  if(task_index == y){
	    found = true;
	    aug_path = _bg.GetTask(y)->path;
	    break;
	  }
	  else{
	    dq.push_back(*_bg.GetTask(task_index));
	    //break;
  	  }
        } //if(!_bg.bg_..)
    }// if "AGENT"

    else{
      //wierd situation
      cerr<<"Error: Some vertices are found not initialized with Obj...Stopped."<<endl;
      exit(-1);
    }

  }//end while

  return aug_path;
}


vector<EID> 
Hungarian::BFSAugmentingPath(VID x, VID y){

  return this->BFSAugmentingPath(this->bg, x, y);

}


bool 
Hungarian::NeedReLabel(vector<VID>& _T, vector<VID>& _N){

  //Below has been done by RefreshBG();
  //get the neighbor _N for _S
  //_N.clear();
  //_N = this->FindNeighbors(_EG, _S);
  
  //then compare the elements in _N and _T
  if(_N.size() != _T.size())
    return false; 
  else{
    sort(_N.begin(), _N.end());
    sort(_T.begin(), _T.end());
    if(_N == _T) return true;
    else return false;
  }

}

bool 
Hungarian::NeedReLabel(void){

  return this->NeedReLabel(this->T, this->N);

}

vector<VID>
Hungarian::FindNeighbors(const vector<EID>& _EG, const vector<VID>& _S){
  vector<VID> _N;
  _N.clear();
  for(vector<EID>::const_iterator e_itr = _EG.begin(); e_itr != _EG.end(); e_itr++)
    for(vector<VID>::const_iterator v1_itr=_S.begin(); v1_itr != _S.end(); v1_itr++)
      if((*e_itr).first == *v1_itr){
	vector<VID>::iterator v2_itr;
        for(v2_itr=_N.begin(); v2_itr!= _N.end(); v2_itr++)
	  if((*e_itr).second == *v2_itr) break;
	if(v2_itr==_N.end())
	  _N.push_back((*e_itr).second);
      } //if
  return _N;
}


vector<VID> 
Hungarian::FindNeighbors(void){

  return this->FindNeighbors(this->EG, this->S);

}



VID 
Hungarian::PickFreeAgent(BipartiteGraph& _bg){

  //if still not perfect
  int free_agent = -1;

  for(size_t i=0; i<_bg.GetNumAgents(); i++)
    if(!_bg.GetAgent(i)->GetMatched()){
      free_agent = i;
      break;
    }

  if(-1 == free_agent){
    cerr<<"Error: No free agent vertex available. Stopped."<<endl;
    exit(-1);
  }

  return (VID)free_agent;
}


VID 
Hungarian::PickFreeAgent(void){
 
  return this->PickFreeAgent(this->bg);

}


VID 
Hungarian::PickAvailableTask(vector<VID>& _T, vector<VID>& _N){
  int y = -1;
  for(vector<VID>::iterator itr1 = _N.begin(); itr1 != _N.end(); itr1++){
    vector<VID>::iterator itr2;
    for(itr2 = _T.begin(); itr2 != _T.end(); itr2++)
      if(*itr2 == *itr1) break;
    if(itr2 == _T.end()){
      y=*itr1;
      break;
    }
  }

  if(-1 == y){
    cerr<<"Error: No free vertex can be picked! Stopped."<<endl; 
    exit(-1);
  }

  return (VID)y; 

}


VID 
Hungarian::PickAvailableTask(void){

  return this->PickAvailableTask(this->T, this->N);

}


double
Hungarian::UpdateLabels(BipartiteGraph& _bg){
  double delta;
  double _min_delta = POS_INF;
  
  //calculate the mininal delta for re-labelling
  for(size_t i=0; i<_bg.GetNumAgents(); i++){
    if(_bg.GetAgent(i)->GetColored())                //agent is in set of S
      for(size_t j=0; j<_bg.GetNumTasks(); j++){
	if(!_bg.GetTask(j)->GetColored()){            //task is NOT in set of T
	  delta = _bg.GetAgent(i)->GetLabel() + _bg.GetTask(j)->GetLabel() - _bg.GetMatrix(i, j)->GetWeight(); 
	  if(delta < _min_delta && fabs(delta - _min_delta) > DOUBLE_EPSILON){ 
	    _min_delta = delta;
	  }
	}//if
      }//interior for
  }//exterior for

  //cout<<"min_delta: "<<_min_delta<<endl;
  _bg.SetMinDelta(_min_delta);

  //agents re-labelling
  for(size_t i=0; i<_bg.GetNumAgents(); i++)
    //agent is in set of S
    if(_bg.GetAgent(i)->GetColored()){
      double lb = _bg.GetAgent(i)->GetLabel() - _min_delta;
      _bg.GetAgent(i)->SetLabel(lb);
    }

  //tasks re-labelling
  for(size_t j=0; j<_bg.GetNumTasks(); j++)
    //task is in set of T
    if(_bg.GetTask(j)->GetColored()){
      double lb = _bg.GetTask(j)->GetLabel() + _min_delta;
      _bg.GetTask(j)->SetLabel(lb);
    }
  
  return _min_delta;
}


double 
Hungarian::UpdateLabels(void){

  return this->UpdateLabels(this->bg);

}


double
Hungarian::OptimalValue(BipartiteGraph& _bg, vector<EID>& _M){

  double sum = 0;
  for(vector<EID>::iterator itr=_M.begin(); itr!=_M.end(); itr++)
    sum += _bg.GetMatrix(*itr)->GetWeight();

  return sum;

}


double
Hungarian::OptimalValue(void){

  return this->OptimalValue(this->bg, this->M);

}



void
Hungarian::DisplayData(vector<VID>& v){
  for(vector<VID>::iterator itr = v.begin(); itr != v.end(); itr++)
    cout<<*itr<<" ";
  cout<<endl;
}

void
Hungarian::DisplayData(vector<EID>& e){
  for(vector<EID>::iterator itr = e.begin(); itr != e.end(); itr++)
    cout<<"("<<(*itr).first<<","<<(*itr).second<<") ";
  cout<<endl;
}

void
Hungarian::DisplayData(bool s, bool t, bool n, bool eg, bool m){

  //displaying S
  if(s){
    cout<<"Set S:"<<endl;
    DisplayData(S); 
  }

  //displaying T
  if(t){
    cout<<"Set T:"<<endl;
    DisplayData(T);
  }

  //displaying N
  if(n){
    cout<<"Set N:"<<endl;
    DisplayData(N);
  }

  //displaying EG
  if(eg){
    cout<<"Set EG:"<<endl;
    DisplayData(EG);
  }

  //displaying M
  if(m){
    cout<<"Set M:"<<endl;
    DisplayData(M);
  }

}


void
Hungarian::DisplayLabels(vector<Vertex>& vertices){
  for(vector<Vertex>::iterator itr = vertices.begin(); itr != vertices.end(); itr++)
    cout<<(*itr).GetLabel()<<" ";
  cout<<endl;
}


void
Hungarian::DisplayMatrix(Matrix& _m){

  if(_m[0].size() > 30){
    cout<<"Matrix is big, no displaying. "<<endl;
    return;
  }

  for(unsigned int i=0; i<_m.size(); i++){
    for(unsigned int j=0; j<_m[0].size(); j++)
      cout<<"  "<<_m[i][j].GetWeight()<<"\t";
    cout<<endl;
  }
  cout<<endl;

}


void
Hungarian::DisplayConfig(BipartiteGraph& _bg){

  if(_bg.GetNumTasks() > 30){
    cout<<"Matrix is big, no displaying. "<<endl;
    return;
  }

  for(unsigned int i=0; i<_bg.GetNumAgents(); i++){
    for(unsigned int j=0; j<_bg.GetNumTasks(); j++){
      if(_bg.GetMatrix(i,j)->GetMatchedFlag())
        cout<<"  1\t";
      else
        cout<<"  0\t";
    }
    cout<<endl;
  }
  cout<<endl;

}




