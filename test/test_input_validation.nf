nextflow.enable.dsl = 2

include { validateSampleSheetRows; validateHisat2IndexFiles } from '../lib/input_validation'

def expectFailure(String expectedMessage, Closure action) {
    try {
        action.call()
    } catch (Throwable error) {
        def cause = error.cause ?: error
        assert cause.message?.contains(expectedMessage): cause.toString()
        return
    }
    throw new IllegalStateException("Expected validation failure containing '${expectedMessage}'.")
}

workflow {
    def inputRoot = file(launchDir)
    def rows = [
        [sample_id: 'test_pos1', fastq_1: "${inputRoot}/test/test_pos1_1.fastq.gz", fastq_2: "${inputRoot}/test/test_pos1_2.fastq.gz", condition: 'positive'],
        [sample_id: 'test_pos2', fastq_1: "${inputRoot}/test/test_pos2_1.fastq.gz", fastq_2: "${inputRoot}/test/test_pos2_2.fastq.gz", condition: 'positive'],
        [sample_id: 'test_neg1', fastq_1: "${inputRoot}/test/test_neg1_1.fastq.gz", fastq_2: "${inputRoot}/test/test_neg1_2.fastq.gz", condition: 'negative'],
        [sample_id: 'test_neg2', fastq_1: "${inputRoot}/test/test_neg2_1.fastq.gz", fastq_2: "${inputRoot}/test/test_neg2_2.fastq.gz", condition: 'negative']
    ]

    assert validateSampleSheetRows(rows, 'negative').size() == 4
    expectFailure('unique') { validateSampleSheetRows(rows + rows[0], 'negative') }
    expectFailure('Exactly two condition labels') { validateSampleSheetRows(rows + [[sample_id: 'other', fastq_1: rows[0].fastq_1, fastq_2: rows[0].fastq_2, condition: 'other']], 'negative') }
    expectFailure('not present') { validateSampleSheetRows(rows, 'missing') }
    expectFailure('does not exist') {
        def missingReadRows = rows.collect { row -> row.clone() }
        missingReadRows[0].fastq_1 = "${inputRoot}/test/missing.fastq.gz"
        validateSampleSheetRows(missingReadRows, 'negative')
    }
    expectFailure('missing required columns') {
        validateSampleSheetRows([[sample_id: 'sample']], 'negative')
    }
    expectFailure('only letters') {
        def unsafeIdRows = rows.collect { row -> row.clone() }
        unsafeIdRows[0].sample_id = 'unsafe;name'
        validateSampleSheetRows(unsafeIdRows, 'negative')
    }

    def smallIndex = (1..8).collect { part -> file("${inputRoot}/test/genome.${part}.ht2") }
    def largeIndex = (1..8).collect { part -> file("${inputRoot}/test/genome.${part}.ht2l") }
    assert validateHisat2IndexFiles(smallIndex, 'genome') == smallIndex
    assert validateHisat2IndexFiles(largeIndex, 'genome') == largeIndex
    expectFailure('incomplete') { validateHisat2IndexFiles(smallIndex.take(7), 'genome') }
    expectFailure('No HISAT2 index files') { validateHisat2IndexFiles([], 'genome') }

    println 'Passed input validation checks: sample sheet fields, IDs, read files, conditions, and both HISAT2 index formats.'
}
