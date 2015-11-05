#ifndef SEQAN_VCF_IO_TABIX_H_
#define SEQAN_VCF_IO_TABIX_H_

#include <cstdlib>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <cstring>
#include <vector>

namespace seqan {

class Tabix
{
 public:
  htsFile* fn;
  tbx_t* tbx;
  hts_itr_t* hts_iter;
  const tbx_conf_t *idxconf;
  int tid = -1;
  int begin = -1;
  int end = -1;
  bool has_jumped = false;
  String<String<char> > chroms;
  int rID = 0;
  // Iterator<String<String<char> > >::Type current_chrom;
  // std::vector<std::string>::iterator current_chrom;
  // std::string filename;
  

  // Tabix(void) { }

  // ~Tabix(void)
  // {
  //   tbx_itr_destroy(hts_iter);
  //   tbx_destroy(tbx);
  // }
};

inline void
clear(Tabix index)
{
    // tbx_itr_destroy(index.hts_iter);
    // tbx_destroy(index.tbx);
    index.tid = -1;
    index.begin = -1;
    index.end = -1;
    index.has_jumped = false;
    clear(index.chroms);
    // clear(index.current_chrom);
}

inline void
getHeader(String<char> & header, Tabix & index)
{
  // clear(header);
  kstring_t str = {0,0,0};

  while ( hts_getline(index.fn, KS_SEP_LINE, &str) >= 0 )
  {
    if ( !str.l || str.s[0] != index.tbx->conf.meta_char )
    {
      break;
    }
    else
    {
      append(header, str.s);
      append(header, "\n");
    }
  }

  // set back to start
  // index.rID = 0;

  // if (index.hts_iter)
  // {
  //   tbx_itr_destroy(index.hts_iter);
  // }

  // index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[rID]));
}

inline bool
setRegion(Tabix & index, const char * region)
{
  tbx_itr_destroy(index.hts_iter);
  index.hts_iter = tbx_itr_querys(index.tbx, region);
  index.has_jumped = true;
  return true;
}

inline bool
atEnd(Tabix & index)
{
  return index.rID == length(index.chroms) - 1;
}

inline bool
getNextLine(String<char> & line, Tabix & index)
{
  kstring_t str = {0,0,0};

  if (index.hts_iter && tbx_itr_next(index.fn, index.tbx, index.hts_iter, &str) >= 0)
  {
    line = str.s;
    return true;
  }
  else if (!atEnd(index))
  {
    // The current rID is finished, move to the next one
    ++index.rID;
    tbx_itr_destroy(index.hts_iter);
    index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[index.rID]));

    if (index.hts_iter && tbx_itr_next(index.fn, index.tbx, index.hts_iter, &str) >= 0)
    {
      line = str.s;
      return true;
    }
  }

  return false;
}

inline bool open(Tabix & index, char const * vcfFilename)
{
  clear(index);
  index.has_jumped = false;
  struct stat stat_tbi,stat_vcf;
  char *fnidx = (char*) calloc(strlen(vcfFilename) + 5, 1);
  strcat(strcpy(fnidx, vcfFilename), ".tbi");

  if (bgzf_is_bgzf(vcfFilename)!=1 )
  {
    SEQAN_FAIL("File '%s' was not identified as bgzipped. Please use bgzip to compress the file.", vcfFilename);
    std::free(fnidx);
  }

  // Common source of errors: new VCF is used with an old index
  stat(fnidx, &stat_tbi);
  stat(vcfFilename, &stat_vcf);

  if (stat_vcf.st_mtime > stat_tbi.st_mtime )
  {
    SEQAN_FAIL("The index file is older than the bcf file. Please reindex the bcf file.");
  }

  std::free(fnidx);

  if ((index.fn = hts_open(vcfFilename, "r")) == 0)
  {
    SEQAN_FAIL("Fail to open the VCF file.");
  }

  if ((index.tbx = tbx_index_load(vcfFilename)) == NULL)
  {
    SEQAN_FAIL("Failed to load the VCF index file.");
  }

  int nseq;
  const char** seq = tbx_seqnames(index.tbx, &nseq);

  for (int i = 0; i < nseq; ++i)
  {
    appendValue(index.chroms, seq[i]);
  }

  std::free(seq);
  index.idxconf = &tbx_conf_vcf;

  // set up the iterator, defaults to the beginning
  index.rID = 0;
  index.hts_iter = tbx_itr_querys(index.tbx, toCString(index.chroms[0]));
}

}  // namespace seqan 

#endif  // SEQAN_VCF_IO_TABIX_H_
