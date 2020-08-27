#include <cstdio>
#include <cstring>

#include <algorithm>
#include <chrono>
#include <numeric>
#include <thread>
#include <unordered_set>
#include <vector>

using namespace std;

static const bool DEF_SHOW_PROGRESS = false;
static const bool DEF_WAIT_KEY = false;
static const int DEF_SLEEP_MS = 100;

static chrono::milliseconds sleeps = chrono::milliseconds(DEF_SLEEP_MS);

typedef unsigned long long bitm;

struct Segment {
    bitm mask;
    int len;
    int min_shift;
    int max_shift;
};

class Nono {
public:
    explicit Nono(vector<vector<int>>& rows, vector<vector<int>>& cols);
    bool solve();
    void show();
    void setOption(bool showProgress, bool waitKey);

private:
    bool show_progress = DEF_SHOW_PROGRESS;
    bool wait_key = DEF_WAIT_KEY;

    bitm pos_masks[64];

    int nrow, ncol;
    bitm full_row, full_col;

    vector<vector<Segment>> segments_row;
    vector<vector<Segment>> segments_col;

    vector<bitm> omasks_row;
    vector<bitm> omasks_col;
    vector<bitm> xmasks_row;
    vector<bitm> xmasks_col;

    bitm common_omask, common_xmask;

    int line_runs;

    void prepareLine(vector<Segment>& dst, vector<int>& src, int limit, int* sum);
    bitm lenToBitMask(int len);

    void markOverlaps();
    void markO(int row, int col);
    void markX(int row, int col);
    bool runLine(vector<Segment>& segments, bitm omask, bitm xmask, int limit);
    bool moveSegment(vector<Segment>& segments, bitm omask, bitm xmask, int idx, int shift_start, int limit, bitm covered, bitm uncovered);
    bitm updateResult(bitm result, int idx, vector<bitm>& lines, vector<bitm>& crosses);

    bool isRowFinished(int row);
    bool isColFinished(int col);

    void showProgress(int row, int col);
    void showInternal(int row, int col);
    char getSymbol(int row, int col);

    static char O, X, U;
};

char Nono::O = '@';
char Nono::X = '=';
char Nono::U = '.';

Nono::Nono(vector<vector<int>>& rows, vector<vector<int>>& cols)
{
    for (auto i = 0; i < 64; i++) {
        pos_masks[i] = static_cast<bitm>(1) << i;
    }

    nrow = static_cast<int>(rows.size());
    ncol = static_cast<int>(cols.size());
    printf("rows: %d, cols: %d\n", nrow, ncol);
    if (nrow > 64 || ncol > 64) {
        throw exception();
    }
    full_row = lenToBitMask(ncol);
    full_col = lenToBitMask(nrow);

    int sum_row = 0, sum_col = 0;
    segments_row = vector<vector<Segment>>(nrow);
    segments_col = vector<vector<Segment>>(ncol);
    for (auto i = 0; i < nrow; i++) {
        prepareLine(segments_row[i], rows[i], ncol, &sum_row);
    }
    for (auto i = 0; i < ncol; i++) {
        prepareLine(segments_col[i], cols[i], nrow, &sum_col);
    }
    if (sum_row != sum_col) {
        printf("Sum of row values %d and col values %d are different\n", sum_row, sum_col);
    }

    omasks_row = vector<bitm>(nrow, 0);
    omasks_col = vector<bitm>(ncol, 0);
    xmasks_row = vector<bitm>(nrow, 0);
    xmasks_col = vector<bitm>(ncol, 0);
}

void Nono::prepareLine(vector<Segment>& dst, vector<int>& src, int limit, int* sum)
{
    auto cnt = src.size();
    if (cnt == 1 && src[0] == 0) {
        return;
    }
    int position = -1;
    for (auto i = 0; i < cnt; i++) {
        position++;
        int len = src[i];
        Segment segment;
        segment.len = len;
        segment.mask = lenToBitMask(len);
        segment.min_shift = position;
        dst.push_back(segment);
        position += len;
        *sum += len;
    }
    int margin = limit - position;
    if (margin < 0) {
        printf("A row with length sum %d does not fit in the column length %d\n", position, limit);
    }
    for (auto i = 0; i < cnt; i++) {
        dst[i].max_shift = dst[i].min_shift + margin;
    }
}

bitm Nono::lenToBitMask(int len)
{
    bitm mask = 0;
    for (auto i = 0; i < len; i++) {
        mask = (mask << 1) | 1;
    }
    return mask;
}

bool Nono::solve()
{
    line_runs = 0;

    markOverlaps();

    bitm changed_row = accumulate(omasks_col.begin(), omasks_col.end(), static_cast<bitm>(0),
        [](bitm x, bitm y) { return x | y; });
    changed_row = accumulate(xmasks_col.begin(), xmasks_col.end(), changed_row,
        [](bitm x, bitm y) { return x | y; });

    bitm changed_col = accumulate(omasks_row.begin(), omasks_row.end(), static_cast<bitm>(0),
        [](bitm x, bitm y) { return x | y; });
    changed_col = accumulate(xmasks_row.begin(), xmasks_row.end(), changed_col,
        [](bitm x, bitm y) { return x | y; });

    bool finished;
    do {
        finished = true;
        for (auto row = 0; row < nrow; row++) {
            if (!(changed_row & pos_masks[row])) {
                finished = finished && isRowFinished(row);
                continue;
            }
            common_omask = full_row;
            common_xmask = full_row;
            if (!runLine(segments_row[row], omasks_row[row], xmasks_row[row], ncol)) {
                printf("Cannot solve this problem 1\n");
                return false;
            }
            auto changed = updateResult(common_omask, row, omasks_row, omasks_col)
                | updateResult(common_xmask, row, xmasks_row, xmasks_col);
            if (changed != 0) {
                changed_col |= changed;
                showProgress(row, -1);
            }
            finished = finished && isRowFinished(row);
            changed_row ^= pos_masks[row];
        }
        for (auto col = 0; col < ncol; col++) {
            if (!(changed_col & pos_masks[col])) {
                finished = finished && isColFinished(col);
                continue;
            }
            common_omask = full_col;
            common_xmask = full_col;
            if (!runLine(segments_col[col], omasks_col[col], xmasks_col[col], nrow)) {
                printf("Cannot solve this problem 2\n");
                return false;
            }
            auto changed = updateResult(common_omask, col, omasks_col, omasks_row)
                | updateResult(common_xmask, col, xmasks_col, xmasks_row);
            if (changed != 0) {
                changed_row |= changed;
                showProgress(-1, col);
            }
            finished = finished && isColFinished(col);
            changed_col ^= pos_masks[col];
        }
        if (!changed_row && !changed_col && !finished) {
            printf("No changed line left\n");
            return false;
        }
    } while (!finished);

    printf("Total line runs: %d\n", line_runs);
    return true;
}

void Nono::markOverlaps()
{
    for (auto row = 0; row < nrow; row++) {
        bool updated = false;
        for (auto& segment : segments_row[row]) {
            auto end = segment.min_shift + segment.len;
            for (auto col = segment.max_shift; col < end; col++) {
                markO(row, col);
                updated = true;
            }
        }
        if (segments_row[row].empty()) {
            for (auto col = 0; col < ncol; col++) {
                markX(row, col);
            }
            updated = true;
        }
        if (updated) showProgress(row, -1);
    }
    for (auto col = 0; col < ncol; col++) {
        bool updated = false;
        for (auto& segment : segments_col[col]) {
            auto end = segment.min_shift + segment.len;
            for (auto row = segment.max_shift; row < end; row++) {
                markO(row, col);
                updated = true;
            }
        }
        if (segments_col[col].empty()) {
            for (auto row = 0; row < nrow; row++) {
                markX(row, col);
            }
            updated = true;
        }
        if (updated) showProgress(-1, col);
    }
}

void Nono::markO(int row, int col)
{
    omasks_row[row] |= pos_masks[col];
    omasks_col[col] |= pos_masks[row];
}

void Nono::markX(int row, int col)
{
    xmasks_row[row] |= pos_masks[col];
    xmasks_col[col] |= pos_masks[row];
}

bool Nono::runLine(vector<Segment>& segments, bitm omask, bitm xmask, int limit)
{
    line_runs++;
    return moveSegment(segments, omask, xmask, 0, 0, limit, 0, 0);
}

bool Nono::moveSegment(vector<Segment>& segments, bitm omask, bitm xmask, int idx, int shift_start, int limit, bitm covered, bitm uncovered)
{
    if (idx == segments.size()) {
        for (auto i = max(shift_start - 1, 0); i < limit; i++) {
            uncovered |= pos_masks[i];
        }
        if (uncovered & omask) {
            return false;
        }
        common_omask &= covered;
        common_xmask &= uncovered;
        return true;
    }
    auto res = false;
    auto& segment = segments[idx];
    auto seg_mask = segment.mask;
    auto seg_len = segment.len;
    for (auto i = shift_start; i <= segment.max_shift; i++) {
        if (i > 0) {
            uncovered |= pos_masks[i - 1];
        }
        if (uncovered & omask) {
            return res;
        }
        bitm new_covered = covered | (seg_mask << i);
        if (new_covered & xmask) {
            continue;
        }
        res |= moveSegment(segments, omask, xmask, idx + 1, i + seg_len + 1, limit, new_covered, uncovered);
    }
    return res;
}

bitm Nono::updateResult(bitm result, int idx, vector<bitm>& lines, vector<bitm>& crosses)
{
    bitm org = lines[idx];
    if (result == org) {
        return 0;
    }
    if ((result | org) != result) {
        printf("Shouldn't be here 1\n");
        throw exception();
    }
    lines[idx] = result;

    bitm changed = result ^ org;
    bitm crossUpdated = pos_masks[idx];
    int limit = static_cast<int>(crosses.size());
    if (limit != nrow && limit != ncol) {
        printf("Shouldn't be here 2\n");
        throw exception();
    }
    for (auto i = 0; i < limit; i++) {
        if (pos_masks[i] & changed) {
            crosses[i] |= crossUpdated;
        }
    }
    return changed;
}

bool Nono::isRowFinished(int row)
{
    return (omasks_row[row] | xmasks_row[row]) == full_row;
}

bool Nono::isColFinished(int col)
{
    return (omasks_col[col] | xmasks_col[col]) == full_col;
}

void Nono::show()
{
    showInternal(-1, -1);
}

void Nono::showProgress(int row, int col)
{
    if (show_progress) {
        showInternal(row, col);
    }
}

void Nono::showInternal(int row, int col)
{
    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" \
        "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    printf("===");
    if (row >= 0) {
        printf(" row %d ", row + 1);
    }
    if (col >= 0) {
        printf(" col %d ", col + 1);
    }
    printf("===\n");
    for (auto col = 0; col < ncol + 2; col++) {
        if (col % 5 == 0) {
            printf("+ ");
        } else {
            printf("- ");
        }
    }
    printf("\n");
    for (auto r = 0; r < nrow; r++) {
        if (r % 5 == 4) {
            printf("+ ");
        } else {
            printf("| ");
        }
        for (auto c = 0; c < ncol; c++) {
            printf("%c ", getSymbol(r, c));
        }
        if (r % 5 == 4) {
            printf("+");
        } else {
            printf("|");
        }
        if (wait_key && col >= 0) {
            printf("  ");
            printf("%c ", getSymbol(r, col));
        }
        printf("\n");
    }
    for (auto col = 0; col < ncol + 2; col++) {
        if (col % 5 == 0) {
            printf("+ ");
        } else {
            printf("- ");
        }
    }
    printf("\n");

    if (wait_key && row >= 0) {
        printf("\n  ");
        for (auto c = 0; c < ncol; c++) {
            printf("%c ", getSymbol(row, c));
        }
        printf("\n");
    } else {
        printf("\n\n");
    }

    if (wait_key) {
        char bb[16];
        fgets(bb, 16, stdin);
    } else {
        this_thread::sleep_for(sleeps);
    }
}

char Nono::getSymbol(int row, int col)
{
    auto col_mask = pos_masks[col];
    auto x = omasks_row[row] & col_mask ? O : U;
    if (x == U) {
        x = xmasks_row[row] & col_mask ? X : U;
    }
    return x;
}

void Nono::setOption(bool showProgress, bool waitKey)
{
    show_progress = showProgress;
    wait_key = waitKey;
}

///////////////////////////////////////////////////////////////////////////////
// Run and test
///////////////////////////////////////////////////////////////////////////////
static bool showProgress = false;
static bool waitKey = false;
static bool longTest = false;

void runCommon(vector<vector<int>>& rows, vector<vector<int>>& cols)
{
    auto start = chrono::system_clock::now();
    Nono nono(rows, cols);
    nono.setOption(showProgress, waitKey);
    bool success = nono.solve();
    auto end = chrono::system_clock::now();
    nono.show();
    if (success) printf("SUCCESS: ");
    else printf("FAILURE: ");
    printf("took %lld us.\n", chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}

void test1()
{
    vector<vector<int>> rows({
        vector<int>({2, 1}),
        vector<int>({2, 3, 3}),
        vector<int>({5, 2, 4}),
        vector<int>({3, 3, 4}),
        vector<int>({1, 3, 2, 3}),
        vector<int>({5, 3, 1, 3, 2}),
        vector<int>({9, 6, 2, 1, 3}),
        vector<int>({1, 1, 3, 1, 3, 1}),
        vector<int>({9, 3, 2, 2, 3, 2}),
        vector<int>({4, 4, 2, 2, 1, 2, 1}),
        vector<int>({5, 3, 6, 2, 3, 2}),
        vector<int>({1, 2, 6, 2, 1, 2, 2}),
        vector<int>({1, 3, 4, 1, 5, 2, 2}),
        vector<int>({4, 3, 3, 2, 2, 2, 1}),
        vector<int>({4, 2, 5, 3, 2, 2}),
        vector<int>({3, 3, 1, 5, 3, 1, 4}),
        vector<int>({2, 2, 1, 1, 5, 7, 3}),
        vector<int>({1, 1, 4, 5, 2, 4, 3}),
        vector<int>({2, 1, 12, 2, 2, 2}),
        vector<int>({1, 1, 1, 3, 3, 2}),
        vector<int>({2, 2, 3, 3, 3}),
        vector<int>({5, 14, 5}),
        vector<int>({5, 6}),
        vector<int>({1, 1, 13}),
        vector<int>({2, 6, 10}),
        vector<int>({2, 5}),
        vector<int>({3, 12}),
        vector<int>({3, 5}),
        vector<int>({14}),
        vector<int>({10}),
        });
    vector<vector<int>> cols({
        vector<int>({12}),
        vector<int>({5, 4, 1, 2}),
        vector<int>({1, 8, 2}),
        vector<int>({2, 7, 2}),
        vector<int>({1, 2, 1, 1, 1, 7}),
        vector<int>({4, 2, 2, 3}),
        vector<int>({1, 2, 8, 2, 3}),
        vector<int>({2, 13, 2}),
        vector<int>({1, 6, 2, 2, 2}),
        vector<int>({2, 3, 2, 4, 2, 2, 2}),
        vector<int>({5, 1, 1, 1, 1, 3}),
        vector<int>({1, 4, 2, 5, 1, 1, 2}),
        vector<int>({5, 3, 7, 1, 1, 2}),
        vector<int>({1, 2, 3, 8, 1, 1, 1, 2}),
        vector<int>({1, 1, 1, 7, 1, 2, 1, 2}),
        vector<int>({1, 2, 5, 1, 1, 1, 2}),
        vector<int>({1, 2, 5, 1, 1, 1, 1, 2}),
        vector<int>({6, 7, 1, 2, 1, 2}),
        vector<int>({2, 1, 4, 1, 2, 1, 2}),
        vector<int>({5, 1, 2, 1, 2}),
        vector<int>({4, 1, 1, 1, 2, 4}),
        vector<int>({1, 3, 5, 2, 7}),
        vector<int>({2, 3, 3, 2, 1, 7}),
        vector<int>({1, 2, 4, 3, 2, 7}),
        vector<int>({3, 1, 4, 2, 3, 7}),
        vector<int>({2, 6, 5, 2, 5}),
        vector<int>({4, 2, 2, 3, 5}),
        vector<int>({7, 2, 7}),
        vector<int>({2, 3, 3, 5}),
        vector<int>({3}),
        });
    runCommon(rows, cols);
}

void test2()
{
    vector<vector<int>> rows({
        vector<int>({2, 2}),
        vector<int>({2, 2}),
        vector<int>({2, 2}),
        vector<int>({2, 2}),
        vector<int>({8}),
        vector<int>({10}),
        vector<int>({10}),
        vector<int>({2, 4, 2}),
        vector<int>({4, 4}),
        vector<int>({8}),
        });
    vector<vector<int>> cols({
        vector<int>({4}),
        vector<int>({6}),
        vector<int>({7, 2}),
        vector<int>({10}),
        vector<int>({4, 1}),
        vector<int>({4, 1}),
        vector<int>({10}),
        vector<int>({7, 2}),
        vector<int>({6}),
        vector<int>({4}),
        });
    runCommon(rows, cols);
}

void test()
{
    if (longTest)
        test1();
    else
        test2();
}

vector<vector<int>> buildLines(FILE* fp, int nlines)
{
    char buf[256];
    vector<vector<int>> result;
    for (int i = 0; i < nlines; i++) {
        if (fgets(buf, 256, fp) == NULL) {
            printf("Failed to read line from file\n");
            return result;
        }
        vector<int> line;
        const char* sep = " \t\n\r(){}[]<>`~!@#$%^&*-_=+\\|;:'\",./?abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
        for (char* tok = strtok(buf, sep); tok; tok = strtok(NULL, sep)) {
            int n;
            sscanf(tok, "%d", &n);
            line.push_back(n);
        }
        printf("%2d > Read %d numbers...", i + 1, static_cast<int>(line.size()));
        for_each(line.begin(), line.end(), [](auto x) { printf(" %d", x); });
        printf("\n");
        result.emplace_back(line);
    }
    return result;
}

void runFile(const char* filename)
{
    char buf[16];
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Cannot open file %s\n", filename);
        return;
    }
    int nrow, ncol;
    if (fgets(buf, 16, fp) == NULL) {
        printf("Failed to read line from %s\n", filename);
        fclose(fp);
        return;
    }
    sscanf(buf, "%d %d", &nrow, &ncol);
    printf("%d rows and %d columns\n", nrow, ncol);
    auto rows = buildLines(fp, nrow);
    if (rows.size() != nrow) {
        fclose(fp);
        return;
    }
    auto cols = buildLines(fp, ncol);
    if (cols.size() != ncol) {
        fclose(fp);
        return;
    }
    fclose(fp);

    runCommon(rows, cols);
}

void setOpt(int argc, const char* argv[])
{
    for (auto i = 1; i < argc; i++) {
        auto a = argv[i];
        if (a[0] != '-')
            continue;
        auto l = strlen(a);
        for (auto j = 1; j < l; j++) {
            auto c = a[j];
            if (c == 's') {
                showProgress = true;
            } else if (c == 'w') {
                waitKey = true;
            } else if (c == 'l') {
                longTest = true;
            }
        }
    }
}

const char* getFilename(int argc, const char* argv[])
{
    for (auto i = 1; i < argc; i++) {
        auto a = argv[i];
        if (a[0] != '-')
            return a;
    }
    return nullptr;
}

int main(int argc, const char* argv[])
{
    setOpt(argc, argv);
    auto name = getFilename(argc, argv);
    if (name == nullptr) {
        test();
    } else {
        runFile(name);
    }
    return 0;
}
