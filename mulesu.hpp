#include <Eigen/Dense>
#include <map>
#include <vector>
#include <regex>
#include <algorithm>

namespace mulesu
{
    namespace const_variables
    {
        std::map<std::string, float> atomicMass = 
        {
            {"H", 1.007975},
            {"He", 4.002602},
            {"Li", 6.967499999999999},
            {"Be", 9.0121831},
            {"B", 10.8135},
            {"C", 12.0106},
            {"N", 14.006855},
            {"O", 15.9994},
            {"F", 18.998403163},
            {"Ne", 20.1797},
            {"Na", 22.98976928},
            {"Mg", 24.3055},
            {"Al", 26.9815385},
            {"Si", 28.085},
            {"P", 30.973761998},
            {"S", 32.067499999999995},
            {"Cl", 35.451499999999996},
            {"Ar", 39.948},
            {"K", 39.0983},
            {"Ca", 40.078},
            {"Sc", 44.955908},
            {"Ti", 47.867},
            {"V", 50.9415},
            {"Cr", 51.9961},
            {"Mn", 54.938044},
            {"Fe", 55.845},
            {"Co", 58.933194},
            {"Ni", 58.6934},
            {"Cu", 63.546},
            {"Zn", 65.38},
            {"Ga", 69.723},
            {"Ge", 72.63},
            {"As", 74.921595},
            {"Se", 78.971},
            {"Br", 79.904},
            {"Kr", 83.798},
            {"Rb", 85.4678},
            {"Sr", 87.62},
            {"Y", 88.90584},
            {"Zr", 91.224},
            {"Nb", 92.90637},
            {"Mo", 95.95},
            {"Tc", 98.0},
            {"Ru", 101.07},
            {"Rh", 102.9055},
            {"Pd", 106.42},
            {"Ag", 107.8682},
            {"Cd", 112.414},
            {"In", 114.818},
            {"Sn", 118.71},
            {"Sb", 121.76},
            {"Te", 127.6},
            {"I", 126.90447},
            {"Xe", 131.293},
            {"Cs", 132.90545196},
            {"Ba", 137.327},
            {"La", 138.90547},
            {"Ce", 140.116},
            {"Pr", 140.90766},
            {"Nd", 144.242},
            {"Pm", 145.0},
            {"Sm", 150.36},
            {"Eu", 151.964},
            {"Gd", 157.25},
            {"Tb", 158.92535},
            {"Dy", 162.5},
            {"Ho", 164.93033},
            {"Er", 167.259},
            {"Tm", 168.93422},
            {"Yb", 173.054},
            {"Lu", 174.9668},
            {"Hf", 178.49},
            {"Ta", 180.94788},
            {"W", 183.84},
            {"Re", 186.207},
            {"Os", 190.23},
            {"Ir", 192.217},
            {"Pt", 195.084},
            {"Au", 196.966569},
            {"Hg", 200.592},
            {"Tl", 204.3835},
            {"Pb", 207.2},
            {"Bi", 208.9804},
            {"Po", 209.0},
            {"At", 210.0},
            {"Rn", 222.0},
            {"Fr", 223.0},
            {"Ra", 226.0},
            {"Ac", 227.0},
            {"Th", 232.0377},
            {"Pa", 231.03588},
            {"U", 238.02891},
            {"Np", 237.0},
            {"Pu", 244.0},
            {"Am", 243.0},
            {"Cm", 247.0},
            {"Bk", 247.0},
            {"Cf", 251.0},
            {"Es", 252.0},
            {"Fm", 257.0},
            {"Md", 258.0},
            {"No", 259.0},
            {"Lr", 262.0},
            {"Rf", 267.0},
            {"Db", 268.0},
            {"Sg", 271.0},
            {"Bh", 272.0},
            {"Hs", 270.0},
            {"Mt", 276.0},
            {"Ds", 281.0},
            {"Rg", 280.0},
            {"Cn", 285.0},
            {"Nh", 284.0},
            {"Fl", 289.0},
            {"Mc", 288.0},
            {"Lv", 293.0},
            {"Ts", 292.0},
            {"Og", 295.0}
        };
        std::regex pat_big("~([0-9|.]{0,})([A-Z|a-z|0-9]{1,})");
        std::regex pat_sml("([A-Z][a-z]?)(\\d{0,})");
    }

    namespace _internal
    {
        class _RegexIter
        {
            public:
                _RegexIter(const std::string& str, const std::regex& pat):
                    _pat(pat), _begin(str.begin(), str.end(), _pat), _end()
                {};
                auto begin()
                {
                    return _begin;
                }
                auto end()
                {
                    return _end;
                }
            private:
                std::regex _pat;
                std::sregex_iterator _begin;
                std::sregex_iterator _end;
        };

        auto calMw(const std::string& eqn)
        {
            float mw = 0;
            std::string effi_str;
            std::string elem;
            int effi;
            for (auto m : _RegexIter(eqn, const_variables::pat_sml))
            {
                effi_str = m[2];
                elem = m[1];
                if (effi_str == "") effi = 1;
                else effi = std::stoi(effi_str);
                mw += effi * const_variables::atomicMass[elem];
            }
            return mw;
        }

        auto _combination(int n, int r)
        {
            assert(n >= r && r > 0);
            std::vector<std::vector<int>> output;
            std::vector<int> buff;
            if (r == 1)
            {
                for (int i = 0; i < n; ++i)
                {
                    buff.clear();
                    buff.push_back(i);
                    output.push_back(buff);
                }
                return output;
            }
            else if (n == r)
            {
                for (int i = 0; i < n; ++i) buff.push_back(i);
                output.push_back(buff);
                return output;
            }
            else
            {
                auto befList = _combination(n-1, r);
                for (auto& v : befList)
                {
                    output.push_back(v);
                }
                befList = _combination(n-1, r-1);
                for (auto& v : befList)
                {
                    buff.clear();
                    buff = v;
                    buff.push_back(n-1);
                    output.push_back(buff);
                }
                return output;
            }
        }

        template<typename rowType, typename colType, typename Scalar>
        class FrameBase
        // 확장 가능한 배열의 기본 형태
        {
            public:
                typedef typename Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic> DynTypeArr;

            protected:
                // FrameBase의 기본적인 멤버 변수들
                std::vector<rowType> _rowStr;
                std::vector<colType> _colStr;
                std::map<rowType, int> _rowMap;
                std::map<colType, int> _colMap;
                DynTypeArr _dataArray;

            public:
                // 생성자 정의부
                FrameBase() = default;
                FrameBase(const DynTypeArr& dataArray, const std::vector<rowType>& rowStr,
                    const std::vector<colType>& colStr):
                    _rowStr(rowStr), _colStr(colStr), _dataArray(dataArray)
                {
                    assert(_rowStr.size() == _dataArray.rows() && _colStr.size() == _dataArray.cols());
                    for (auto i = 0; i < _rowStr.size(); ++i) _rowMap[_rowStr[i]] = i;
                    for (auto j = 0; j < _colStr.size(); ++j) _colMap[_colStr[j]] = j;
                };
                FrameBase(const FrameBase& other):
                    _rowStr(other._rowStr), _colStr(other._colStr),
                    _dataArray(other._dataArray)
                {
                    assert(_rowStr.size() == _dataArray.rows() && _colStr.size() == _dataArray.cols());
                    for (auto i = 0; i < _rowStr.size(); ++i) _rowMap[_rowStr[i]] = i;
                    for (auto j = 0; j < _colStr.size(); ++j) _colMap[_colStr[j]] = j;
                };

                FrameBase& operator+=(const FrameBase& other)
                {
                    std::vector<rowType> r_rowStr = _rowStr;
                    std::vector<colType> r_colStr = _colStr;
                    DynTypeArr r_dataArray;

                    typename std::map<rowType, int>::const_iterator s_rowIter;
                    typename std::map<rowType, int>::const_iterator o_rowIter;
                    for (auto i = 0; i < other._rowStr.size(); ++i)
                    {
                        s_rowIter = _rowMap.find(other._rowStr[i]);
                        if (s_rowIter == _rowMap.end()) r_rowStr.push_back(other._rowStr[i]);
                    }

                    typename std::map<colType, int>::const_iterator s_colIter;
                    typename std::map<colType, int>::const_iterator o_colIter;
                    for (auto j = 0; j < other._colStr.size(); ++j)
                    {
                        s_colIter = _colMap.find(other._colStr[j]);
                        if (s_colIter == _colMap.end()) r_colStr.push_back(other._colStr[j]);
                    }

                    r_dataArray.resize(r_rowStr.size(), r_colStr.size());
                    r_dataArray.setZero();
                    
                    for (auto i = 0; i < r_rowStr.size(); ++i)
                    {   
                        s_rowIter = _rowMap.find(r_rowStr[i]);
                        o_rowIter = other._rowMap.find(r_rowStr[i]);

                        for (auto j = 0; j < r_colStr.size(); ++j)
                        {
                            s_colIter = _colMap.find(r_colStr[j]);
                            o_colIter = other._colMap.find(r_colStr[j]);
                            
                            if (s_rowIter != _rowMap.end() && s_colIter != _colMap.end())
                            {
                                r_dataArray(i, j) += _dataArray(s_rowIter->second, s_colIter->second);
                            }

                            if (o_rowIter != other._rowMap.end() && o_colIter != other._colMap.end())
                            {
                                r_dataArray(i, j) += other._dataArray(o_rowIter->second, o_colIter->second);
                            }
                        }
                    }

                    _dataArray = r_dataArray;
                    _rowStr = r_rowStr;
                    _colStr = r_colStr;
                    for (auto i = 0; i < _rowStr.size(); ++i) _rowMap[_rowStr[i]] = i;
                    for (auto j = 0; j < _colStr.size(); ++j) _colMap[_colMap[j]] = j;

                    return *this;
                };
                FrameBase& operator-=(const FrameBase& other)
                {
                    std::vector<rowType> r_rowStr = _rowStr;
                    std::vector<colType> r_colStr = _colStr;
                    DynTypeArr r_dataArray;

                    typename std::map<rowType, int>::const_iterator s_rowIter;
                    typename std::map<rowType, int>::const_iterator o_rowIter;
                    for (auto i = 0; i < other._rowStr.size(); ++i)
                    {
                        s_rowIter = _rowMap.find(other._rowStr[i]);
                        if (s_rowIter == _rowMap.end()) r_rowStr.push_back(other._rowStr[i]);
                    }

                    typename std::map<colType, int>::const_iterator s_colIter;
                    typename std::map<colType, int>::const_iterator o_colIter;
                    for (auto j = 0; j < other._colStr.size(); ++j)
                    {
                        s_colIter = _colMap.find(other._colStr[j]);
                        if (s_colIter == _colMap.end()) r_colStr.push_back(other._colStr[j]);
                    }

                    r_dataArray.resize(r_rowStr.size(), r_colStr.size());
                    r_dataArray.setZero();
                    
                    for (auto i = 0; i < r_rowStr.size(); ++i)
                    {   
                        s_rowIter = _rowMap.find(r_rowStr[i]);
                        o_rowIter = other._rowMap.find(r_rowStr[i]);

                        for (auto j = 0; j < r_colStr.size(); ++j)
                        {
                            s_colIter = _colMap.find(r_colStr[j]);
                            o_colIter = other._colMap.find(r_colStr[j]);
                            
                            if (s_rowIter != _rowMap.end() && s_colIter != _colMap.end())
                            {
                                r_dataArray(i, j) += _dataArray(s_rowIter->second, s_colIter->second);
                            }

                            if (o_rowIter != other._rowMap.end() && o_colIter != other._colMap.end())
                            {
                                r_dataArray(i, j) -= other._dataArray(o_rowIter->second, o_colIter->second);
                            }
                        }
                    }

                    _dataArray = r_dataArray;
                    _rowStr = r_rowStr;
                    _colStr = r_colStr;
                    for (auto i = 0; i < _rowStr.size(); ++i) _rowMap[_rowStr[i]] = i;
                    for (auto j = 0; j < _colStr.size(); ++j) _colMap[_colMap[j]] = j;

                    return *this;
                };
                FrameBase& operator*=(const float& rhs)
                {
                    _dataArray *= rhs;
                    return *this;
                }
                FrameBase& operator/=(const float& rhs)
                {
                    _dataArray /= rhs;
                    return *this;
                }
                FrameBase operator+(const FrameBase& other)
                {
                    FrameBase r_FrameBase(*this);
                    r_FrameBase += other;
                    return r_FrameBase;
                };
                FrameBase operator-(const FrameBase& other)
                {
                    FrameBase r_FrameBase(*this);
                    r_FrameBase -= other;
                    return r_FrameBase;
                };
                FrameBase operator*(const float& rhs)
                {
                    FrameBase r_FrameBase(*this);
                    r_FrameBase *= rhs;
                    return r_FrameBase;
                };
                FrameBase operator/(const float& rhs)
                {
                    FrameBase r_FrameBase(*this);
                    r_FrameBase /= rhs;
                    return r_FrameBase;
                };

                // FrameBase의 getter
                auto getRowStr() {return _rowStr;}
                auto getColStr() {return _colStr;}
                auto getDataArray() {return _dataArray;}

                // 인덱싱을 위한 클래스(모든 행 또는 열을 의미)
                struct All{};

                // 입력한 행 또는 열 이름의 인덱스를 반환함(유효하지 않은 경우 런타임 에러 발생)
                std::vector<int> _getRowIdx(const std::vector<rowType>& rowStr)
                {
                    typename std::map<rowType, int>::const_iterator rowIter;
                    std::vector<int> rowIdx;
                    for (auto& rowName : rowStr)
                    {
                        rowIter = _rowMap.find(rowName);
                        if (rowIter == _rowMap.end()) throw std::runtime_error("invalid row index");
                        rowIdx.push_back(rowIter->second);
                    }
                    return rowIdx;
                }
                int _getRowIdx(const rowType& rowName)
                {
                    typename std::map<rowType, int>::const_iterator rowIter = _rowMap.find(rowName);
                    if (rowIter == _rowMap.end()) throw std::runtime_error("invalid row index");
                    return rowIter->second;
                }
                std::vector<int> _getColIdx(const std::vector<colType>& colStr)
                {
                    typename std::map<colType, int>::const_iterator colIter;
                    std::vector<int> colIdx;
                    for (auto& colName : colStr)
                    {
                        colIter = _colMap.find(colName);
                        if (colIter == _colMap.end()) throw std::runtime_error("invalid col index");
                        colIdx.push_back(colIter->second);
                    }
                    return colIdx;
                }
                int _getColIdx(const colType& colName)
                {
                    typename std::map<colType, int>::const_iterator colIter = _colMap.find(colName);
                    if (colIter == _colMap.end()) throw std::runtime_error("invalid col index");
                    return colIter->second;
                }
                
                // FrameBase의 인덱싱
                Scalar operator()(const rowType& rowStr, const colType& colStr)
                {
                    auto rowIdx = _getRowIdx(rowStr);
                    auto colIdx = _getColIdx(colStr);
                    return _dataArray(rowIdx, colIdx);
                }
                FrameBase operator()(const std::vector<rowType>& rowStr, const All&)
                {
                    auto rowIdx = _getRowIdx(rowStr);
                    DynTypeArr r_dataArray(rowStr.size(), _colStr.size());
                    for (auto i = 0; i < r_dataArray.rows(); ++i)
                    {
                        r_dataArray.row(i) = _dataArray.row(rowIdx[i]);
                    }
                    FrameBase r_FrameBase(r_dataArray, rowStr, _colStr);
                    return r_FrameBase;
                }
                FrameBase operator()(const All&, const std::vector<colType>& colStr)
                {
                    auto colIdx = _getColIdx(colStr);
                    DynTypeArr r_dataArray(_rowStr.size(), colstr.size());
                    for (auto j = 0; j < r_dataArray.cols(); ++j)
                    {
                        r_dataArray.col(j) = _dataArray.row(colIdx[j]);
                    }
                    FrameBase r_FrameBase(r_dataArray, _rowStr, colStr);
                    return r_FrameBase;
                }
                FrameBase operator()(const std::vector<rowType>& rowStr, const std::vector<colType>& colStr)
                {
                    auto rowIdx = _getRowIdx(rowStr);
                    auto colIdx = _getColIdx(colStr);
                    DynTypeArr r_data
                }
                auto row(const rowType& rowStr)
                {
                    typename std::map<rowType, int>::const_iterator rowIter = _rowMap.find(rowStr);
                    if (rowIter == _rowMap.end())
                    {
                        throw std::runtime_error("row index error");
                    }
                    return _dataArray.row(rowIter->second);
                }

                auto col(const rowType& colStr)
                {
                    typename std::map<colType, int>::const_iterator colIter = _colMap.find(colStr);
                    if (colIter == _colMap.end())
                    {
                        throw std::runtime_error("col index error");
                    }
                    return _dataArray.col(colIter->second);
                }

        };
    }

    
}