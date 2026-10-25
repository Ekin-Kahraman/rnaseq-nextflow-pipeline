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
        [sample_id: 'test_pos1', fastq_1: "${inputRoot}/test/test_pos1_1.fastq.gz", fastq_2: "${inputRoot}/test/test_pos1_2.fastq.gz", condition: 'positive', donor: 'donor_1'],
        [sample_id: 'test_pos2', fastq_1: "${inputRoot}/test/test_pos2_1.fastq.gz", fastq_2: "${inputRoot}/test/test_pos2_2.fastq.gz", condition: 'positive', donor: 'donor_2'],
        [sample_id: 'test_neg1', fastq_1: "${inputRoot}/test/test_neg1_1.fastq.gz", fastq_2: "${inputRoot}/test/test_neg1_2.fastq.gz", condition: 'negative', donor: 'donor_1'],
        [sample_id: 'test_neg2', fastq_1: "${inputRoot}/test/test_neg2_1.fastq.gz", fastq_2: "${inputRoot}/test/test_neg2_2.fastq.gz", condition: 'negative', donor: 'donor_2']
    ]

    assert validateSampleSheetRows(rows, 'negative').size() == 4
    assert validateSampleSheetRows(rows, 'negative', 'donor_condition').size() == 4
    expectFailure('Unsupported DESeq2 design') { validateSampleSheetRows(rows, 'negative', 'arbitrary_formula') }
    expectFailure('no residual degrees of freedom') {
        validateSampleSheetRows([rows[0], rows[2]], 'negative', 'condition')
    }
    assert validateSampleSheetRows(rows.take(3), 'negative', 'condition').size() == 3
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
    expectFailure('leading or trailing whitespace') {
        def paddedIdRows = rows.collect { row -> row.clone() }
        paddedIdRows[0].sample_id = ' test_pos1'
        validateSampleSheetRows(paddedIdRows, 'negative')
    }
    expectFailure('Condition labels must not have leading or trailing whitespace') {
        def paddedConditionRows = rows.collect { row -> row.clone() }
        paddedConditionRows[0].condition = 'positive '
        validateSampleSheetRows(paddedConditionRows, 'negative')
    }
    expectFailure('missing required columns') {
        def missingDonorRows = rows.collect { row -> row.findAll { key, _value -> key != 'donor' } }
        validateSampleSheetRows(missingDonorRows, 'negative', 'donor_condition')
    }
    expectFailure('at least two donors') {
        def singleDonorRows = rows.collect { row ->
            def copy = row.clone()
            copy.donor = 'donor_1'
            copy
        }
        validateSampleSheetRows(singleDonorRows, 'negative', 'donor_condition')
    }
    expectFailure('Every donor must have samples in both conditions') {
        def unpairedRows = rows.collect { row -> row.clone() }
        unpairedRows[0].donor = 'donor_2'
        validateSampleSheetRows(unpairedRows, 'negative', 'donor_condition')
    }
    expectFailure('Donor labels must not have leading or trailing whitespace') {
        def paddedDonorRows = rows.collect { row -> row.clone() }
        paddedDonorRows[0].donor = ' donor_1'
        validateSampleSheetRows(paddedDonorRows, 'negative', 'donor_condition')
    }

    def smallIndex = (1..8).collect { part -> file("${inputRoot}/test/genome.${part}.ht2") }
    def largeIndex = (1..8).collect { part -> file("${inputRoot}/test/genome.${part}.ht2l") }
    assert validateHisat2IndexFiles(smallIndex, 'genome') == smallIndex
    assert validateHisat2IndexFiles(largeIndex, 'genome') == largeIndex
    def neighbours = (1..8).collect { part -> file("${inputRoot}/test/genome-old.${part}.ht2") }
    def largeNeighbours = (1..8).collect { part -> file("${inputRoot}/test/genome_other.${part}.ht2l") }
    assert validateHisat2IndexFiles(neighbours + smallIndex.reverse(), 'genome') == smallIndex
    assert validateHisat2IndexFiles(largeNeighbours + largeIndex.reverse(), 'genome') == largeIndex
    assert validateHisat2IndexFiles(largeIndex + neighbours + smallIndex, 'genome') == smallIndex
    expectFailure('incomplete') { validateHisat2IndexFiles(neighbours + smallIndex.take(7), 'genome') }
    expectFailure('Ambiguous') { validateHisat2IndexFiles(smallIndex + smallIndex[0], 'genome') }
    expectFailure('incomplete') { validateHisat2IndexFiles(smallIndex.take(7), 'genome') }
    expectFailure('No HISAT2 index files') { validateHisat2IndexFiles([], 'genome') }

    def literalRows = rows.collect { row -> row.clone() }
    literalRows.eachWithIndex { row, i ->
        row.sample_id = ['001', '1', 'NA', '002'][i]
        row.donor = ['01', '1', '01', '1'][i]
    }
    assert validateSampleSheetRows(literalRows, 'negative', 'donor_condition') == literalRows

    println 'Passed input validation checks: sample sheet fields, IDs, paired donor design, read files, conditions, and both HISAT2 index formats.'
}
