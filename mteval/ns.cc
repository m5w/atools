#include "ns.h"
#include "ns_ter.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;
using boost::shared_ptr;

map<string, EvaluationMetric*> EvaluationMetric::instances_;

SegmentEvaluator::~SegmentEvaluator() {}
EvaluationMetric::~EvaluationMetric() {}

struct DefaultSegmentEvaluator : public SegmentEvaluator {
  DefaultSegmentEvaluator(const vector<vector<WordID> >& refs, const EvaluationMetric* em) : refs_(refs), em_(em) {}
  void Evaluate(const vector<WordID>& hyp, SufficientStats* out) const {
    em_->ComputeSufficientStatistics(hyp, refs_, out);
  }
  const vector<vector<WordID> > refs_;
  const EvaluationMetric* em_;
};

shared_ptr<SegmentEvaluator> EvaluationMetric::CreateSegmentEvaluator(const vector<vector<WordID> >& refs) const {
  return shared_ptr<SegmentEvaluator>(new DefaultSegmentEvaluator(refs, this));
}

void EvaluationMetric::ComputeSufficientStatistics(const vector<WordID>&,
                                                   const vector<vector<WordID> >&,
                                                   SufficientStats*) const {
  cerr << "Base class ComputeSufficientStatistics should not be called.\n";
  abort();
}

enum BleuType { IBM, Koehn, NIST };
template <unsigned int N = 4u, BleuType BrevityType = IBM>
struct BleuSegmentEvaluator : public SegmentEvaluator {
  BleuSegmentEvaluator(const vector<vector<WordID> >& refs, const EvaluationMetric* em) : evaluation_metric(em) {
    assert(refs.size() > 0);
    float tot = 0;
    int smallest = 9999999;
    for (vector<vector<WordID> >::const_iterator ci = refs.begin();
         ci != refs.end(); ++ci) {
      lengths_.push_back(ci->size());
      tot += lengths_.back();
      if (lengths_.back() < smallest) smallest = lengths_.back();
      CountRef(*ci);
    }
    if (BrevityType == Koehn)
      lengths_[0] = tot / refs.size();
    if (BrevityType == NIST)
      lengths_[0] = smallest;
  }

  void Evaluate(const vector<WordID>& hyp, SufficientStats* out) const {
    out->fields.resize(N + N + 2);
    out->evaluation_metric = evaluation_metric;
    for (unsigned i = 0; i < N+N+2; ++i) out->fields[i] = 0;

    ComputeNgramStats(hyp, &out->fields[0], &out->fields[N], true);
    float& hyp_len = out->fields[2*N];
    float& ref_len = out->fields[2*N + 1];
    hyp_len = hyp.size();
    ref_len = lengths_[0];
    if (lengths_.size() > 1 && BrevityType == IBM) {
      float bestd = 2000000;
      float hl = hyp.size();
      float bl = -1;
      for (vector<float>::const_iterator ci = lengths_.begin(); ci != lengths_.end(); ++ci) {
        if (fabs(*ci - hl) < bestd) {
          bestd = fabs(*ci - hl);
          bl = *ci;
        }
      }
      ref_len = bl;
    }
  }

  struct NGramCompare {
    int operator() (const vector<WordID>& a, const vector<WordID>& b) {
      const size_t as = a.size();
      const size_t bs = b.size();
      const size_t s = (as < bs ? as : bs);
      for (size_t i = 0; i < s; ++i) {
         int d = a[i] - b[i];
         if (d < 0) return true;
         if (d > 0) return false;
      }
      return as < bs;
    }
  };
  typedef map<vector<WordID>, pair<int,int>, NGramCompare> NGramCountMap;

  void CountRef(const vector<WordID>& ref) {
    NGramCountMap tc;
    vector<WordID> ngram(N);
    int s = ref.size();
    for (int j=0; j<s; ++j) {
      int remaining = s-j;
      int k = (N < remaining ? N : remaining);
      ngram.clear();
      for (int i=1; i<=k; ++i) {
        ngram.push_back(ref[j + i - 1]);
        tc[ngram].first++;
      }
    }
    for (typename NGramCountMap::iterator i = tc.begin(); i != tc.end(); ++i) {
      pair<int,int>& p = ngrams_[i->first];
      if (p.first < i->second.first)
        p = i->second;
    }
  }

  void ComputeNgramStats(const vector<WordID>& sent,
                         float* correct,  // N elements reserved
                         float* hyp,      // N elements reserved
                         bool clip_counts = true) const {
    vector<WordID> ngram(N);
    *correct *= 0;
    *hyp *= 0;
    int s = sent.size();
    for (int j=0; j<s; ++j) {
      int remaining = s-j;
      int k = (N < remaining ? N : remaining);
      ngram.clear();
      for (int i=1; i<=k; ++i) {
        ngram.push_back(sent[j + i - 1]);
        pair<int,int>& p = ngrams_[ngram];
        if(clip_counts){
          if (p.second < p.first) {
            ++p.second;
            correct[i-1]++;
          }
        } else {
          ++p.second;
          correct[i-1]++;
        }
        // if the 1 gram isn't found, don't try to match don't need to match any 2- 3- .. grams:
        if (!p.first) {
          for (; i<=k; ++i)
            hyp[i-1]++;
        } else {
          hyp[i-1]++;
        }
      }
    }
  }

  const EvaluationMetric* evaluation_metric;
  vector<float> lengths_;
  mutable NGramCountMap ngrams_;
};

template <unsigned int N = 4u, BleuType BrevityType = IBM>
struct BleuMetric : public EvaluationMetric {
  BleuMetric() : EvaluationMetric("IBM_BLEU") {}
  float ComputeScore(const SufficientStats& stats) const {
    float log_bleu = 0;
    int count = 0;
    for (int i = 0; i < N; ++i) {
      if (stats.fields[i+N] > 0) {
        float cor_count = stats.fields[i];  // correct_ngram_hit_counts[i];
        // smooth bleu
        if (!cor_count) { cor_count = 0.01; }
        float lprec = log(cor_count) - log(stats.fields[i+N]); // log(hyp_ngram_counts[i]);
        // if (precs) precs->push_back(exp(lprec));
        log_bleu += lprec;
        ++count;
      }
    }
    log_bleu /= count;
    float lbp = 0.0;
    const float& hyp_len = stats.fields[2*N];
    const float& ref_len = stats.fields[2*N + 1];
    if (hyp_len < ref_len)
      lbp = (hyp_len - ref_len) / hyp_len;
    log_bleu += lbp;
    //if (bp) *bp = exp(lbp);
    return exp(log_bleu);
  }
  shared_ptr<SegmentEvaluator> CreateSegmentEvaluator(const vector<vector<WordID> >& refs) const {
    return shared_ptr<SegmentEvaluator>(new BleuSegmentEvaluator<N,BrevityType>(refs, this));
  }
};

EvaluationMetric* EvaluationMetric::Instance(const string& metric_id) {
  static bool is_first = true;
  if (is_first) {
    instances_["NULL"] = NULL;
    is_first = false;
  }

  map<string, EvaluationMetric*>::iterator it = instances_.find(metric_id);
  if (it == instances_.end()) {
    EvaluationMetric* m = NULL; 
    if (metric_id == "IBM_BLEU") {
      m = new BleuMetric<4, IBM>;
    } else if (metric_id == "NIST_BLEU") {
      m = new BleuMetric<4, NIST>;
    } else if (metric_id == "Koehn_BLEU") {
      m = new BleuMetric<4, Koehn>;
    } else if (metric_id == "TER") {
      m = new TERMetric;
    } else {
      cerr << "Implement please: " << metric_id << endl;
      abort();
    }
    if (m->MetricId() != metric_id) {
      cerr << "Registry error: " << metric_id << " vs. " << m->MetricId() << endl;
      abort();
    }
    return instances_[metric_id] = m;
  } else {
    return it->second;
  }
}

SufficientStats::SufficientStats(const string& encoded) {
  istringstream is(encoded);
  string type;
  is >> type;
  evaluation_metric = EvaluationMetric::Instance(type);
  float val;
  while(is >> val)
    fields.push_back(val);
}

void SufficientStats::Encode(string* out) const {
  ostringstream os;
  if (evaluation_metric)
    os << evaluation_metric->MetricId();
  else
    os << "NULL";
  for (unsigned i = 0; i < fields.size(); ++i)
    os << ' ' << fields[i];
  *out = os.str();
}

