#include <vector>
#include <map>
#include <fstream>
#include <stdlib.h>
using namespace std;

void split(const string& str, vector<string>, char del);
typedef map<uint32_t, uint32_t> TopicCountDist;

struct Document{
    vector<uint32_t> _words;
    vector<uint32_t> _topics;
    TopicCountDist _doc_count_dist;
};

struct Model{
    vector<TopicCountDist*> _topic_count_dist;
    vector<uint32_t> _topic_count_total;
    uint32_t _topic_num;
};

double rand_double()
{
    return rand() / static_cast<double>(RAND_MAX);
}
uint32_t rand_uint(uint32_t bound)
{
    return static_cast<uint32_t>(rand_double()) * bound;
}

//class Corpus {
//public:
//    Corpus() {}
//    ~Corpus() {}
//
//private:
//    vector<Document> _corpus;
//    map<uint32_t, string> _id2word;
//    map<string, uint32_t> _word2id;
//};
//
class AliasSampler {
public:
    AliasSampler(const Model& model, const vector<Document>& corpus, double alpha, double beta, uint32_t voc_num)
    : _model(model), _corpus(corpus), _alpha(alpha), _beta(beta), _voc_num(voc_num) {}
    ~AliasSampler() {}

    uint32_t do_sample(uint32_t docidx, uint32_t wordidx)
    {
        Document& doc = corpus[docidx];
        vector<pair<uint32_t, double> > doc_topic_dist;
        double p_sum = 0.0;
        for (TopicCountDist::iterator it=doc._doc_count_dist.begin(); it!=doc._doc_count_dist.end(); ++it)
        {
            uint32_t topic = it->first;
            uint32_t word = doc._words[wordidx];

            uint32_t adjusted_doc_topic_num = it->second;
            uint32_t adjusted_topic_word_num = model._topic_count_total[it->first];
            uint32_t adjusted_total_topic_word_num = 0;
            if (model._topic_count_dist[topic]->find(word) != model._topic_count_dist[topic]->end())
                adjusted_total_topic_word_num = (*model._topic_count_dist[topic])[word];
            
            if (topic == doc._topics[wordidx])
            {
                adjusted_doc_topic_num --; 
                adjusted_topic_word_num --;
                adjusted_total_topic_word_num --;
            }
            double prob = adjusted_doc_topic_num * (adjusted_topic_word_num + _beta)   
                          / (adjusted_total_topic_word_num + _voc_num*_beta);

            doc_topic_dist.push_back(make_pair(topic, prob)); 
            p_sum += prob;
        }   
        
        for (uint32_t i=0; i<doc_topic_dist.size(); ++i)
            doc_topic_dist[i].second /= p_sum;   

        
    }
private:
    struct AliasTable {
        uint32_t l;
        uint32_t h;
        double prob;
    };
    
    void build_alias_table(uint32_t topic, vector<AliasTable>& alias_table)
    {
        uint32_t topic_count_total = _model._topic_count_total[topic];

        list<pair<uint32_t, double> > L;
        list<pair<uint32_t, double> > H;       
        double avg = 1.0 / _voc_num;
        for (uint32_t i=0; i<_voc_num; ++i)
        {
            double prob = 0.0; 
            if (TopicCountDist.find(i) != TopicCountDist.end())
                prob = (it->second + beta) / (model._topic_num + beta*_voc_num); 
            else
                prob = (it->second + beta) / (model.topic_num + V*beta); 

            if (prob <= avg)
                L.push_back(make_pair(it->first, prob));
            else
                H.push_back(make_pair(it->first, prob));
        }

        while (L.size()>0 && H.size()>0)
        {
            pair<uint32_t, double> l = L.front(); 
            pair<uint32_t, double> h = H.front(); 
            L.pop_front();
            H.pop_fornt();
            AliasTable at;
            at.l = l.first;
            at.h = h.first;
            at.prob = l.second;
            alias_table.push_back(at);
            h.second = l.second + h.second - avg;
            if (h.second <= avg)
                L.push_back(h);
            else
                H.push_back(h);
        }
       
        list<pair<uint32_t, double> >& left = L.size()!=0 ? L : H;      
        while (left.size() != 0)
        {
            AliasTable at;
            pair<uint32_t, double> p = left.front();
            left.pop_front();
            at.l = p.first;
            at.h = 0;
            at.prob = p.second; 
        } 
    }

    void sample_from_alias_table(uint32_t num, const vector<AliasTable>& alias_table)
    {
        for (uint32_t i=0; i<num; ++i)
        {
            uint32_t t = rand_uint(_voc_num);
            double prob = rand_double();
            if (prob/_voc_num <= alias_table[t].prob)
                _candidate_sample.push_back(alias_table[t].l);
            else
                _candidate_sample.push_back(alias_table[t].h);
        }
    }


private:
    list<uint32_t> _candidate_sample;
    double _alpha;
    double _beta;
    uint32_t _voc_num;
    uint32_t _topic_num; 
    
    Model& _model; 
    vector<Document>& _corpus;
 //   double generate_alias_sample(
};

class LDATrainer {
public:
    LDATrainer(uint32_t topic_num, double alpha, double beta, const string& corpus_path)
    : _corpus_path(corpus_path), _topic_num(topic_num), _alpha(alpha), _beta(beta) {}
    ~LDATrainer() {}

    void init() 
    {
        for (size_t k=0; k<_topic_num; ++k)
            _model._topic_count_dist[k] = new TopicCountDist();
        _model._topic_count_total.resize(_topic_num);  
        // for each document in corpus
        for (size_t i=0; i<_corpus.size(); ++i)
        {
            // for each word in doc
            Document& doc = _corpus[i];
            for (size_t j=0; j<doc._topics.size(); ++j)
            {
                uint32_t topic = rand_uint(_topic_num);
                doc._topics[j] = topic; 
                if (doc._doc_count_dist.find(topic) != doc._doc_count_dist.end())
                    doc._doc_count_dist[topic] = 0;
                doc._doc_count_dist[topic] ++ ;
                
                TopicCountDist* topic_count_dist = _model._topic_count_dist[topic];
                uint32_t word = doc._words[j];
                if (topic_count_dist->find(word) != topic_count_dist->end())
                    (*topic_count_dist)[word] = 0;
                (*topic_count_dist)[word] ++;
                _model._topic_count_total[topic] ++; 
            }
        }
    }

    void train()
    {
        load_corpus();
        init();
        
        int max_step = 100;
        for (int step=0; step<max_step; ++step)
        {
            for (uint32_t d=0; d<_corpus.size(); ++d)
            {
                Document& doc = _corpus[d];
                uint32_t doc_len = doc._words.size();
                for (uint32_t i=0; i<doc_len; ++i)
                {
                    uint32_t topic = do_sample(_alpha, _beta, doc, _model);
                } 
            } 
        }  
    }

    void load_corpus()
    {
        ifstream ifs(_corpus_path.c_str());
        string buff;
        uint32_t word_idx = 0;
        while (getline(ifs, buff))
        {
            Document doc;
            vector<string> tmp;
            split(buff, tmp, ' ');
            for (size_t i=0; i<tmp.size(); ++i)
            {
                string& word = tmp[i];
                uint32_t id = 0;
                if (_word2id.find(word) == _word2id.end())
                {
                    id = word_idx++;
                    _word2id[word] = id;
                    _id2word[id] = word;
                }
                else
                {
                    id = _word2id[word]; 
                }
                doc._words.push_back(id);
            } 
            doc._topics.resize(doc._words.size());
            _corpus.push_back(doc);
        }     
    }


private:
    string _corpus_path;
     
    vector<Document> _corpus;
    map<uint32_t, string> _id2word;
    map<string, uint32_t> _word2id;
    Model _model; 
    uint32_t _topic_num;
    double _alpha;
    double _beta; 
};


int main()
{
}
