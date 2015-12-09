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

class HtsFile
{
  public:
    const char * filename;  /** @brief The filename of the current file. */
    htsFile * fp;           /** @brief Pointer to the file. */
    bam_hdr_t * hdr;        /** @brief The header of the current file. */
    bam1_t * hts_record;    /** @brief The current HTS record. */
    hts_idx_t * hts_index;  /** @brief The index of the file. */
    hts_itr_t * hts_iter;   /** @brief An iterator that iterates through a certain region in the HTS file. */
    const char * file_mode; /** @brief Which file mode to use. E.g. "r" for reading and "wb" for writing binaries. */

    /**
     * @brief Constructs a new HtsFile object.
     * 
     * @param f The filename of the file.
     * @param mode The file mode to use when opening the file.
     * @return A new HtsFile object.
     */
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

    /**
     * @brief Destructs an HtsFile object.
     */
    ~HtsFile()
    {
        bam_hdr_destroy(hdr);
        hts_close(fp);
    }
};

/**
 * @brief Copies a header from a source HTS file and replaces the header of the target.
 * 
 * @param target Target HTS file.
 * @param source Source HTS file.
 */
inline void
copyHeader(HtsFile & target, HtsFile const & source)
{
    target.hdr = bam_hdr_dup(source.hdr);
}

/**
 * @brief Copies a record from a source HTS file and replaces the record of the target.
 * 
 * @param target Target HTS file.
 * @param source Source HTS file.
 */
inline void
copyRecord(HtsFile & target, HtsFile const & source)
{
    target.hts_record = bam_dup1(source.hts_record);
}

/**
 * @brief Writes a HTS header to disk.
 * 
 * @param file HTS file to get the header from.
 * @returns True on success, otherwise false.
 */
inline bool
writeHeader(HtsFile & file)
{
    return !sam_hdr_write(file.fp, file.hdr);
}

/**
 * @brief Writes a HTS record to disk. Returns 
 * 
 * @param file HTS file to get the record from.
 * @returns True on success, otherwise false.
 */
inline bool
writeRecord(HtsFile & file)
{
    return !sam_write1(file.fp, file.hdr, file.hts_record);
}

/**
 * @brief Loads an index for a HTS file using the default filename.
 * 
 * @param file HTS file to load index for.
 * @returns True on success, otherwise false.
 */
inline bool
loadIndex(HtsFile & file)
{
    file.hts_index = sam_index_load(file.fp, file.filename);
    return file.hts_index != NULL;
}

/**
 * @brief Loads an index for a HTS file with a specific filename.
 * 
 * @param file HTS file to load index for.
 * @param indexFileName The filename of the index.
 */
inline bool
loadIndex(HtsFile & file, const char * indexFileName)
{
    file.hts_index = sam_index_load2(file.fp, file.filename, indexFileName);
    return file.hts_index != NULL;
}

/**
 * @brief Builds an index for BAM or CRAM files using the default filename.
 * 
 * @param file The file to build index for.
 * @param min_shift Force a certain minimum amount of shift. (I think) smaller shifts mean more accurate queries at the cost of index size.
 *                  The default value is 0, which means the default value of htslib will be used.
 * @returns True on success, otherwise false.
 */
inline bool
buildIndex(HtsFile & file, int min_shift = 0)
{
    return !sam_index_build(file.filename, min_shift);
}

/**
 * @brief Builds an index for BAM or CRAM files using a specific filename.
 * 
 * @param file The file to build index for.
 * @param min_shift Force a certain minimum amount of shift. (I think) smaller shifts mean more accurate queries at the cost of index size.
 *                  The default value is 0, which means the default value of htslib will be used.
 * @returns True on success, otherwise false.
 */
inline bool
buildIndex(HtsFile & file, const char * indexFileName, int min_shift = 0)
{
    return !sam_index_build2(file.filename, indexFileName, min_shift);
}

/**
 * @brief Uses the index to go to a certain region of the HTS file.
 * 
 * @param file HTS file to change index on.
 * @param region The region to go to. Should be on one of these formats: chrX, chrX:A, or chrX:A-B.
 */
inline void
setRegion(HtsFile & file, const char * region)
{
    if (file.hts_iter != NULL)
        sam_itr_destroy(file.hts_iter);

    file.hts_iter = sam_itr_querys(file.hts_index, file.hdr, region);
}

/**
 * @brief Read the next record from a HTS file.
 * 
 * @param file HTS file to read from.
 * @returns True on success, otherwise false.
 */
inline bool
readRecord(HtsFile & file)
{
    return sam_read1(file.fp, file.hdr, file.hts_record) >= 0;
}

/**
 * @brief Read the next record from a HTS file and parse it to a sequence record.
 * 
 * @param record Sequencing record to write to.
 * @param file HTS file to read from.
 * @returns True on success, otherwise false.
 */
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

/**
 * @brief Read the next record from a region and parse it to a sequence record.
 * 
 * @param record Sequencing record to write to.
 * @param file HTS file to read from.
 * @returns True on success, otherwise false.
 */
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
