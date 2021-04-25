#include <cstdio>
#include <cstring>

#include <algorithm>
#include <chrono>
#include <numeric>
#include <thread>
#include <unordered_set>
#include <vector>

#define USE_HEURISTIC_INIT

using namespace std;

static const bool DEF_SHOW_PROGRESS = false;
static const bool DEF_WAIT_KEY = false;
static const int DEF_SLEEP_MS = 100;

static chrono::milliseconds sleeps = chrono::milliseconds(DEF_SLEEP_MS);

// Up to 64 bit.
typedef unsigned long long BitMask;

// Segment represents contiguous filled points.
// It includes the length, bitmask and movable range
// in a row or a column.
struct Segment {
    BitMask mask;
    int len;
    int min_shift;
    int max_shift;
};

class Nono {
public:
    explicit Nono(vector<vector<int>>& rows, vector<vector<int>>& cols);
    bool Solve();
    void Show();
    void SetOption(bool show_progress, bool wait_key);

private:
    // Options
    bool show_progress_ = DEF_SHOW_PROGRESS;
    bool wait_key_ = DEF_WAIT_KEY;

    // pos_masks_[n] = 1 << n
    BitMask pos_masks_[64];

    int num_row_, num_col_;
    BitMask full_row_, full_col_;

    vector<vector<Segment>> segments_row_;
    vector<vector<Segment>> segments_col_;

    // Updated mask for O and X marks.
    // This is the progress of solving the puzzle.
    // These value only increases bitwise
    // until omask+xmask fills the row or column.
    // The omasks_row_ and omasks_col_ pair or
    // the xmasks_row_ and xmasks_col_ pair have
    // redundant data, so they should be updated
    // at the same function.
    vector<BitMask> omasks_row_;
    vector<BitMask> omasks_col_;
    vector<BitMask> xmasks_row_;
    vector<BitMask> xmasks_col_;

    // Temporary mask for processing a row or a column.
    BitMask common_omask_, common_xmask_;

    int line_runs_;

    void PrepareLine(vector<Segment>& dst, vector<int>& src, int limit, int* sum);
    BitMask LenToBitMask(int len);

#ifdef USE_HEURISTIC_INIT
    void MarkOverlaps();
    void MarkO(int row, int col);
    void MarkX(int row, int col);
#endif
    bool RunLine(vector<Segment>& segments, BitMask omask, BitMask xmask, int limit);
    bool MoveSegment(vector<Segment>& segments, BitMask omask, BitMask xmask, int idx, int shift_start, int limit, BitMask covered, BitMask uncovered);
    BitMask UpdateResult(BitMask result, int idx, vector<BitMask>& lines, vector<BitMask>& crosses);

    bool IsRowFinished(int row);
    bool IsColFinished(int col);

    void ShowProgress(int row, int col);
    void ShowInternal(int row, int col);
    char GetSymbol(int row, int col);

    static char o_, x_, u_;
};

char Nono::o_ = '@';
char Nono::x_ = '=';
char Nono::u_ = '.';

Nono::Nono(vector<vector<int>>& rows, vector<vector<int>>& cols)
{
    for (auto i = 0; i < 64; i++) {
        pos_masks_[i] = static_cast<BitMask>(1) << i;
    }

    num_row_ = static_cast<int>(rows.size());
    num_col_ = static_cast<int>(cols.size());
    printf("rows: %d, cols: %d\n", num_row_, num_col_);
    if (num_row_ > 64 || num_col_ > 64) {
        throw exception();
    }
    full_row_ = LenToBitMask(num_col_);
    full_col_ = LenToBitMask(num_row_);

    int sum_row = 0, sum_col = 0;
    segments_row_ = vector<vector<Segment>>(num_row_);
    segments_col_ = vector<vector<Segment>>(num_col_);
    for (auto i = 0; i < num_row_; i++) {
        PrepareLine(segments_row_[i], rows[i], num_col_, &sum_row);
    }
    for (auto i = 0; i < num_col_; i++) {
        PrepareLine(segments_col_[i], cols[i], num_row_, &sum_col);
    }
    if (sum_row != sum_col) {
        printf("Sum of row values %d and col values %d are different\n", sum_row, sum_col);
    }

    omasks_row_ = vector<BitMask>(num_row_, 0);
    omasks_col_ = vector<BitMask>(num_col_, 0);
    xmasks_row_ = vector<BitMask>(num_row_, 0);
    xmasks_col_ = vector<BitMask>(num_col_, 0);
}

// Build segments from integer values.
void Nono::PrepareLine(vector<Segment>& dst, vector<int>& src, int limit, int* sum)
{
    int cnt = static_cast<int>(src.size());
    if (cnt == 1 && src[0] == 0) {
        return;
    }
    int position = 0;
    for (auto i = 0; i < cnt; i++) {
        int len = src[i];
        Segment segment;
        segment.len = len;
        segment.mask = LenToBitMask(len);
        segment.min_shift = position;
        dst.push_back(segment);
        position += len + 1;  // Including minimum space.
        *sum += len;
    }
    position--;  // Remove the last space.
    int margin = limit - position;
    if (margin < 0) {
        printf("Sum of the segments %d does not fit in the row or column length %d\n", position, limit);
    }
    for (auto i = 0; i < cnt; i++) {
        dst[i].max_shift = dst[i].min_shift + margin;
    }
}

BitMask Nono::LenToBitMask(int len)
{
    BitMask mask = 0;
    for (auto i = 0; i < len; i++) {
        mask = (mask << 1) | 1;
    }
    return mask;
}

bool Nono::Solve()
{
    line_runs_ = 0;
 
#ifdef USE_HEURISTIC_INIT
    MarkOverlaps();

    BitMask changed_row = accumulate(omasks_col_.begin(), omasks_col_.end(), static_cast<BitMask>(0),
        [](BitMask x, BitMask y) { return x | y; });
    changed_row = accumulate(xmasks_col_.begin(), xmasks_col_.end(), changed_row,
        [](BitMask x, BitMask y) { return x | y; });

    BitMask changed_col = accumulate(omasks_row_.begin(), omasks_row_.end(), static_cast<BitMask>(0),
        [](BitMask x, BitMask y) { return x | y; });
    changed_col = accumulate(xmasks_row_.begin(), xmasks_row_.end(), changed_col,
        [](BitMask x, BitMask y) { return x | y; });
#else
    BitMask changed_row = full_row_;
    BitMask changed_col = full_col_;
#endif

    bool finished;
    do {
        finished = true;
        for (auto row = 0; row < num_row_; row++) {
            if (!(changed_row & pos_masks_[row])) {
                // No need to update unchanged row.
                finished = finished && IsRowFinished(row);
                continue;
            }
            common_omask_ = full_row_;  // Updated in RunLine()
            common_xmask_ = full_row_;  // Updated in RunLine()
            if (!RunLine(segments_row_[row], omasks_row_[row], xmasks_row_[row], num_col_)) {
                printf("Cannot solve this problem 1\n");
                return false;
            }
            auto changed = UpdateResult(common_omask_, row, omasks_row_, omasks_col_)
                | UpdateResult(common_xmask_, row, xmasks_row_, xmasks_col_);
            if (changed != 0) {
                changed_col |= changed;
                ShowProgress(row, -1);
            }
            finished = finished && IsRowFinished(row);
            changed_row ^= pos_masks_[row];  // Clear bit for this row.
        }
        for (auto col = 0; col < num_col_; col++) {
            if (!(changed_col & pos_masks_[col])) {
                // No need to update unchanged column.
                finished = finished && IsColFinished(col);
                continue;
            }
            common_omask_ = full_col_;  // Updated in RunLine()
            common_xmask_ = full_col_;  // Updated in RunLine()
            if (!RunLine(segments_col_[col], omasks_col_[col], xmasks_col_[col], num_row_)) {
                printf("Cannot solve this problem 2\n");
                return false;
            }
            auto changed = UpdateResult(common_omask_, col, omasks_col_, omasks_row_)
                | UpdateResult(common_xmask_, col, xmasks_col_, xmasks_row_);
            if (changed != 0) {
                changed_row |= changed;
                ShowProgress(-1, col);
            }
            finished = finished && IsColFinished(col);
            changed_col ^= pos_masks_[col];  // Clear bit for this column.
        }
        if (!changed_row && !changed_col && !finished) {
            printf("No changed line left\n");
            return false;
        }
    } while (!finished);

    printf("Total line runs: %d\n", line_runs_);
    return true;
}

#ifdef USE_HEURISTIC_INIT
// Heuristic initial marking for faster solution.
void Nono::MarkOverlaps()
{
    for (auto row = 0; row < num_row_; row++) {
        bool updated = false;
        for (auto& segment : segments_row_[row]) {
            auto end = segment.min_shift + segment.len;
            for (auto col = segment.max_shift; col < end; col++) {
                MarkO(row, col);
                updated = true;
            }
        }
        if (segments_row_[row].empty()) {
            for (auto col = 0; col < num_col_; col++) {
                MarkX(row, col);
            }
            updated = true;
        }
        if (updated) ShowProgress(row, -1);
    }
    for (auto col = 0; col < num_col_; col++) {
        bool updated = false;
        for (auto& segment : segments_col_[col]) {
            auto end = segment.min_shift + segment.len;
            for (auto row = segment.max_shift; row < end; row++) {
                MarkO(row, col);
                updated = true;
            }
        }
        if (segments_col_[col].empty()) {
            for (auto row = 0; row < num_row_; row++) {
                MarkX(row, col);
            }
            updated = true;
        }
        if (updated) ShowProgress(-1, col);
    }
}

void Nono::MarkO(int row, int col)
{
    omasks_row_[row] |= pos_masks_[col];
    omasks_col_[col] |= pos_masks_[row];
}

void Nono::MarkX(int row, int col)
{
    xmasks_row_[row] |= pos_masks_[col];
    xmasks_col_[col] |= pos_masks_[row];
}
#endif  // USE_HEURISTIC_INIT

bool Nono::RunLine(vector<Segment>& segments, BitMask omask, BitMask xmask, int limit)
{
    line_runs_++;
    if (segments.empty()) {
        common_omask_ = 0;
        return true;
    }
    return MoveSegment(segments, omask, xmask, 0, segments[0].min_shift, limit, 0, 0);
}

// Move segment to all the possible positions.
// Update common_omask_ and common_xmask_.
// Recursion. O(limit^segments.size())
bool Nono::MoveSegment(vector<Segment>& segments, BitMask omask, BitMask xmask, int idx, int shift_start, int limit, BitMask covered, BitMask uncovered)
{
    if (idx == segments.size()) {
        for (auto i = max(shift_start - 1, 0); i < limit; i++) {
            uncovered |= pos_masks_[i];
        }
        if (uncovered & omask) {
            return false;
        }
        common_omask_ &= covered;
        common_xmask_ &= uncovered;
        return true;
    }
    auto res = false;
    auto& segment = segments[idx];
    auto seg_mask = segment.mask;
    auto seg_len = segment.len;
    for (auto i = shift_start; i <= segment.max_shift; i++) {
        if (i > 0) {
            uncovered |= pos_masks_[i - 1];
        }
        if (uncovered & omask) {
            return res;
        }
        BitMask new_covered = covered | (seg_mask << i);
        if (new_covered & xmask) {
            continue;
        }
        res |= MoveSegment(segments, omask, xmask, idx + 1, i + seg_len + 1, limit, new_covered, uncovered);
    }
    return res;
}

// Update the lines and crosses with result.
// Returns the bitmask of changed points.
BitMask Nono::UpdateResult(BitMask result, int idx, vector<BitMask>& lines, vector<BitMask>& crosses)
{
    BitMask org = lines[idx];
    if (result == org) {
        return 0;
    }
    if ((result | org) != result) {
        printf("Shouldn't be here 1\n");
        throw exception();
    }
    lines[idx] = result;

    BitMask changed = result ^ org;
    BitMask cross_updated = pos_masks_[idx];
    int limit = static_cast<int>(crosses.size());
    if (limit != num_row_ && limit != num_col_) {
        printf("Shouldn't be here 2\n");
        throw exception();
    }
    for (auto i = 0; i < limit; i++) {
        if (pos_masks_[i] & changed) {
            crosses[i] |= cross_updated;
        }
    }
    return changed;
}

bool Nono::IsRowFinished(int row)
{
    return (omasks_row_[row] | xmasks_row_[row]) == full_row_;
}

bool Nono::IsColFinished(int col)
{
    return (omasks_col_[col] | xmasks_col_[col]) == full_col_;
}

void Nono::Show()
{
    ShowInternal(-1, -1);
}

void Nono::ShowProgress(int row, int col)
{
    if (show_progress_) {
        ShowInternal(row, col);
    }
}

void Nono::ShowInternal(int row, int col)
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
    for (auto c = 0; c < num_col_ + 2; c++) {
        if (c % 5 == 0) {
            printf("+ ");
        } else {
            printf("- ");
        }
    }
    printf("\n");
    for (auto r = 0; r < num_row_; r++) {
        if (r % 5 == 4) {
            printf("+ ");
        } else {
            printf("| ");
        }
        for (auto c = 0; c < num_col_; c++) {
            printf("%c ", GetSymbol(r, c));
        }
        if (r % 5 == 4) {
            printf("+");
        } else {
            printf("|");
        }
        if (wait_key_ && col >= 0) {
            printf("  ");
            printf("%c ", GetSymbol(r, col));
        }
        printf("\n");
    }
    for (auto c = 0; c < num_col_ + 2; c++) {
        if (c % 5 == 0) {
            printf("+ ");
        } else {
            printf("- ");
        }
    }
    printf("\n");

    if (wait_key_ && row >= 0) {
        printf("\n  ");
        for (auto c = 0; c < num_col_; c++) {
            printf("%c ", GetSymbol(row, c));
        }
        printf("\n");
    } else {
        printf("\n\n");
    }

    if (wait_key_) {
        char bb[16];
        fgets(bb, 16, stdin);
    } else {
        this_thread::sleep_for(sleeps);
    }
}

char Nono::GetSymbol(int row, int col)
{
    auto col_mask = pos_masks_[col];
    auto x = omasks_row_[row] & col_mask ? o_ : u_;
    if (x == u_) {
        x = xmasks_row_[row] & col_mask ? x_ : u_;
    }
    return x;
}

void Nono::SetOption(bool show_progress, bool wait_key)
{
    show_progress_ = show_progress;
    wait_key_ = wait_key;
}

///////////////////////////////////////////////////////////////////////////////
// Run and test
///////////////////////////////////////////////////////////////////////////////
static bool opt_show_progress = false;
static bool opt_wait_key = false;
static bool opt_long_sample = false;

void RunCommon(vector<vector<int>>& rows, vector<vector<int>>& cols)
{
    auto start = chrono::system_clock::now();
    Nono nono(rows, cols);
    nono.SetOption(opt_show_progress, opt_wait_key);
    bool success = nono.Solve();
    auto end = chrono::system_clock::now();
    nono.Show();
    if (success) printf("SUCCESS: ");
    else printf("FAILURE: ");
    printf("took %lld us.\n", chrono::duration_cast<std::chrono::microseconds>(end - start).count());
}

void RunLongSample()
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
    RunCommon(rows, cols);
}

void RunShortSample()
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
    RunCommon(rows, cols);
}

void RunSample()
{
    if (opt_long_sample)
        RunLongSample();
    else
        RunShortSample();
}

vector<vector<int>> BuildLines(FILE* fp, int nlines)
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

void RunFile(const char* filename)
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
    auto rows = BuildLines(fp, nrow);
    if (rows.size() != nrow) {
        fclose(fp);
        return;
    }
    auto cols = BuildLines(fp, ncol);
    if (cols.size() != ncol) {
        fclose(fp);
        return;
    }
    fclose(fp);

    RunCommon(rows, cols);
}

void SetOpt(int argc, const char* argv[])
{
    for (auto i = 1; i < argc; i++) {
        auto a = argv[i];
        if (a[0] != '-')
            continue;
        auto l = strlen(a);
        for (auto j = 1; j < l; j++) {
            auto c = a[j];
            if (c == 's') {
                opt_show_progress = true;
            } else if (c == 'w') {
                opt_wait_key = true;
            } else if (c == 'l') {
                opt_long_sample = true;
            }
        }
    }
}

const char* GetFilename(int argc, const char* argv[])
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
    SetOpt(argc, argv);
    auto name = GetFilename(argc, argv);
    if (name == nullptr) {
        RunSample();
    } else {
        RunFile(name);
    }
    return 0;
}
