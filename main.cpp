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
    uint32_t do_sample(double alpha, double beta, const Document& doc, const Model& model)
    {

    }
private:
    struct AliasTable {
        uint32_t l;
        uint32_t h;
        double prob;
    };
    void build_alias_table(const TopicCountDist& topic_count_dist, vector<AliasTable>& alias_table)
    {
      
    }
    vector<uint32_t> candidate_samle;
    double q_sum;
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
