#ifndef SEQAN_HTS_IO_HTS_FILE_IN_H_
#define SEQAN_HTS_IO_HTS_FILE_IN_H_

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
// #include <htslib/vcf.h>
#include <htslib/bgzf.h>

#include <seqan/hts_io/hts_alignment_record.h>


namespace seqan
{

struct Hts_;
typedef Tag<Hts_> Hts;
typedef FormattedFile<Hts, Input> HtsFileTest;

class HtsFile
{
  public:
    const char * filename;
    htsFile * fp;
    bam_hdr_t * hdr;
    bam1_t * hts_record;
    hts_idx_t * hts_index;
    hts_itr_t * hts_iter;
    const char * file_mode;

    HtsFile(const char * f, const char * mode = "r")
      : filename(f), fp(NULL), hdr(NULL), hts_record(NULL), hts_index(NULL), hts_iter(NULL), file_mode(mode)
    {
        static const char * read_mode = "r";
        fp = hts_open(filename, file_mode);

        if (fp == NULL)
            SEQAN_FAIL("Could not open file with filename %s", filename);

        if (file_mode == read_mode)
            hdr = sam_hdr_read(fp);

        hts_record = bam_init1();
    }

    ~HtsFile()
    {
        bam_hdr_destroy(hdr);
        hts_close(fp);
    }
};

inline void
copyHeader(HtsFile & target, HtsFile const & source)
{
    target.hdr = bam_hdr_dup(source.hdr);
}

inline void
copyRecord(HtsFile & target, HtsFile const & source)
{
    target.hts_record = bam_dup1(source.hts_record);
}

inline bool
writeHeader(HtsFile & file)
{
    return !sam_hdr_write(file.fp, file.hdr);
}

inline bool
writeRecord(HtsFile & file)
{
    return !sam_write1(file.fp, file.hdr, file.hts_record);
}

inline bool
loadIndex(HtsFile & file)
{
    file.hts_index = sam_index_load(file.fp, file.filename);
    return file.hts_index != NULL;
}

/**
 * @brief Loads an index with a specific filename.
 * 
 * @param file [description]
 * @param indexFileName [description]
 */
inline bool
loadIndex(HtsFile & file, const char * indexFileName)
{
    file.hts_index = sam_index_load2(file.fp, file.filename, indexFileName);
    return file.hts_index != NULL;
}


/**
 * @brief Builds an index for BAM or CRAM files.
 * 
 * @param file The file to build index for.
 * @param min_shift Force a certain minimum amount of shift. (I think) smaller shifts mean more accurate queries at the cost of index size.
 *                  The default value is 0, which means the default value of htslib will be used.
 */
inline bool
buildIndex(HtsFile & file, int min_shift = 0)
{
    return !sam_index_build(file.filename, min_shift);
}

inline bool
buildIndex(HtsFile & file, const char * indexFileName, int min_shift = 0)
{
    return !sam_index_build2(file.filename, indexFileName, min_shift);
}

inline void
setRegion(HtsFile & file, const char * region)
{
    if (file.hts_iter != NULL)
        sam_itr_destroy(file.hts_iter);

    file.hts_iter = sam_itr_querys(file.hts_index, file.hdr, region);
}

inline bool
readRecord(HtsFile & file)
{
    return sam_read1(file.fp, file.hdr, file.hts_record) >= 0;
}

inline bool
readRecord(HtsSequenceRecord & record, HtsFile & file)
{
    if (readRecord(file))
    {
        record.parse(file.hts_record);
        return true;
    }
    else
    {
        // We've reached the end of the file, or an error occured.
        return false;
    }
}

inline bool
readRegion(HtsSequenceRecord & record, HtsFile & file)
{
    if (sam_itr_next(file.fp, file.hts_iter, file.hts_record) >= 0)
    {
        record.parse(file.hts_record);
        return true;
    }
    else
    {
        // We've reached the end of the file, or an error occured.
        return false;
    }
}

} // namespace seqan


#endif  // SEQAN_HTS_IO_HTS_FILE_IN_H_
