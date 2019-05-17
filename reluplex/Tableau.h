/*********************                                                        */
/*! \file Tableau.h
 ** \verbatim
 ** Top contributors (to current version):
 **   Guy Katz
 ** This file is part of the Reluplex project.
 ** Copyright (c) 2016-2017 by the authors listed in the file AUTHORS
 ** (in the top-level source directory) and their institutional affiliations.
 ** All rights reserved. See the file COPYING in the top-level source
 ** directory for licensing information.\endverbatim
 **/

#ifndef __Tableau_h__
#define __Tableau_h__

#include "Map.h"
#include "Debug.h"
#include "FloatUtils.h"
#include "Vector.h"

class Tableau
{
public:
    class Entry
    {
    public:
        Entry()
            : _nextInRow( NULL )
            , _prevInRow( NULL )
            , _nextInColumn( NULL )
            , _prevInColumn( NULL )
            , _value( 0.0 )
        {
        }

        Entry *nextInRow()
        {
            return _nextInRow;
        }

        Entry *prevInRow()
        {
            return _prevInRow;
        }

        Entry *nextInColumn()
        {
            return _nextInColumn;
        }

        Entry *prevInColumn()
        {
            return _prevInColumn;
        }

        const Entry *nextInRow() const
        {
            return _nextInRow;
        }

        const Entry *prevInRow() const
        {
            return _prevInRow;
        }

        const Entry *nextInColumn() const
        {
            return _nextInColumn;
        }

        const Entry *prevInColumn() const
        {
            return _prevInColumn;
        }

        void setNextInRow( Entry *entry )
        {
            _nextInRow = entry;
        }

        void setPrevInRow( Entry *entry )
        {
            _prevInRow = entry;
        }

        void setNextInColumn( Entry *entry )
        {
            _nextInColumn = entry;
        }

        void setPrevInColumn( Entry *entry )
        {
            _prevInColumn = entry;
        }

        void setRow( unsigned row )
        {
            _row = row;
        }

        void setColumn( unsigned column )
        {
            _column = column;
        }

        unsigned getRow() const
        {
            return _row;
        }

        unsigned getColumn() const
        {
            return _column;
        }

        void setValue( double value )
        {
            _value = value;
        }

        double getValue() const
        {
            return _value;
        }

        Entry *_nextInRow;
        Entry *_prevInRow;
        Entry *_nextInColumn;
        Entry *_prevInColumn;

        unsigned _row;
        unsigned _column;

        double _value;
    };

    Tableau( unsigned size )
        : _size( size )
        , _rows( NULL )
        , _columns( NULL )
    {
        _rows = new Tableau::Entry *[size];
        _columns = new Tableau::Entry *[size];

        for ( unsigned i = 0; i < _size; ++i )
        {
            _rows[i] = NULL;
            _columns[i] = NULL;

            _rowSize.append( 0 );
            _columnSize.append( 0 );
        }
    }

    unsigned getNumVars() const
    {
        return _size;
    }

    ~Tableau()
    {
        deleteAllEntries();
        delete []_rows;
        delete []_columns;
    }

    unsigned totalSize() const
    {
        unsigned total = 0;

        for ( unsigned i = 0; i < _size; ++i )
            total += _rowSize.get( i );
        return total;
    }

    void deleteAllEntries()
    {
        for ( unsigned i = 0; i < _size; ++i )
        {
            Entry *current, *next;
            current = _rows[i];

            while ( current )
            {
                next = current->nextInRow();
                delete current;
                current = next;
            }
        }

        for ( unsigned i = 0; i < _size; ++i )
        {
            _rows[i] = NULL;
            _columns[i] = NULL;
            _columnSize[i] = 0;
            _rowSize[i] = 0;
        }
    }

    // By row
    double getCell( unsigned row, unsigned column ) const
    {
        Entry *rowEntry = _rows[row];

        while ( rowEntry != NULL )
        {
            if ( rowEntry->getColumn() == column )
                return rowEntry->getValue();

            rowEntry = rowEntry->nextInRow();
        }

        return 0.0;
    }

    void addEntry( unsigned row, unsigned column, const double &value )
    {
        // We assume that the entry does not currently exist (we don't try to unlink it)
        if ( FloatUtils::isZero( value ) )
            return;

		// 设置这个实体本身的值
        Entry *entry = new Entry;
        entry->setRow( row );
        entry->setColumn( column );
        entry->setValue( value );

		// 设置与它有关系的值，
		// setNextInRow需要的形参是Entry指针，而_rows是指针数组，初始化时已经分配好了空间和地址
		// 那么_rows[row]就是表示一个Entry指针
		// rows[x]初始化时指向null，所以第一个Entry的nextRow和nextColumn都是null

		// 总结：（1）到（2）之间的内容，就是手工生成双向链表（采用头加法），用于判断Tableau中的某一row和column是否还有下一项数据
		// (1)
        entry->setNextInRow( _rows[row] );
        entry->setNextInColumn( _columns[column] );

		// 如果next的Entry不为空，则同样设置关联信息，将其的prev指向当前新生成的Entry
		// 实际上这就是手工双向链表的操作，如果前面的next为null，则表示指向结束

        if ( _rows[row] != NULL )
            _rows[row]->setPrevInRow( entry );
        if ( _columns[column] != NULL )
            _columns[column]->setPrevInColumn( entry );

		//将当前新生成的Entry放入相应的row和column指针数组中，作用是留个记录，
		// 等待下次有同样row或column值的Entry生成时，可以将其添加到prev和next的指针指向中
        _rows[row] = entry;
        _columns[column] = entry;
		// (2)
		

		// 记录下当前行和列有多少个有效值
        ++_rowSize[row];
        ++_columnSize[column];
    }

    bool activeRow( unsigned row ) const
    {
        return _rows[row] != NULL;
    }

    bool activeColumn( unsigned column ) const
    {
        return _columns[column] != NULL;
    }

    void eraseEntry( Entry *entry )
    {
        if ( entry->nextInRow() )
            entry->nextInRow()->setPrevInRow( entry->prevInRow() );
        if ( entry->prevInRow() )
            entry->prevInRow()->setNextInRow( entry->nextInRow() );

        if ( entry->nextInColumn() )
            entry->nextInColumn()->setPrevInColumn( entry->prevInColumn() );
        if ( entry->prevInColumn() )
            entry->prevInColumn()->setNextInColumn( entry->nextInColumn() );

        if ( _rows[entry->getRow()] == entry )
            _rows[entry->getRow()] = entry->nextInRow();
        if ( _columns[entry->getColumn()] == entry )
            _columns[entry->getColumn()] = entry->nextInColumn();

        --_rowSize[entry->getRow()];
        --_columnSize[entry->getColumn()];

        delete entry;
    }

    void eraseRow( unsigned row )
    {
        Entry *entry = _rows[row];
        Entry *next;

        while ( entry != NULL )
        {
            if ( entry->nextInColumn() )
                entry->nextInColumn()->setPrevInColumn( entry->prevInColumn() );
            if ( entry->prevInColumn() )
                entry->prevInColumn()->setNextInColumn( entry->nextInColumn() );
            if ( _columns[entry->getColumn()] == entry )
                _columns[entry->getColumn()] = entry->nextInColumn();

            --_columnSize[entry->getColumn()];

            next = entry->nextInRow();
            delete entry;
            entry = next;
        }

        _rows[row] = NULL;
        _rowSize[row] = 0;
    }

    void eraseColumn( unsigned column )
    {
        Entry *entry = _columns[column];
        Entry *next;

        while ( entry != NULL )
        {
            if ( entry->nextInRow() )
                entry->nextInRow()->setPrevInRow( entry->prevInRow() );
            if ( entry->prevInRow() )
                entry->prevInRow()->setNextInRow( entry->nextInRow() );
            if ( _rows[entry->getRow()] == entry )
                _rows[entry->getRow()] = entry->nextInRow();

            --_rowSize[entry->getRow()];

            next = entry->nextInColumn();
            delete entry;
            entry = next;
        }

        _columns[column] = NULL;
        _columnSize[column] = 0;
    }

	// source对应basic,target对应non-basic
    void addScaledRow( unsigned source, double scale, unsigned target,
                       // We usually want to guarantee that a certain entry has a certain value
                       unsigned guaranteeIndex, double guaranteeValue,
                       unsigned *numCalcs = NULL )
    {
        if ( !activeRow( source ) )
            return;

        _denseMap.clear();

		// 遍历target（pivot中的non-basic，或是pivot之后要更改的其他行）行中出现的每一个Entry，存入denseMap，用于后续判断
		// 如果有存在，就不用生成新Entry,只需要改变值即可，否则就需要生成

		// 个人理解：在pivot的non-basic调用时，是肯定不会有已存在的值的，只有在pivot之后更改其他行时，这里才会有值加入denseMap
        Entry *targetEntry = _rows[target];
        while ( targetEntry != NULL )
        {
            _denseMap[targetEntry->getColumn()] = targetEntry;
            targetEntry = targetEntry->nextInRow();
        }
		
		// 遍历source（原basic或是pivot之后修改其他值时传入的non-basic）中每一行的系数，开始更改值
        Entry *sourceEntry = _rows[source];
        Entry *current;
        unsigned column;
        while ( sourceEntry != NULL )
        {
			// 如果是处理新basic行，current表示遍历原basic行中存在的每一个Entry
			// 如果是处理pivot之后与新basic相关的其他行，current表示遍历新basic行中存在的每一个Entry
            current = sourceEntry;
            sourceEntry = sourceEntry->nextInRow();

            column = current->getColumn();
            Entry *entryInTarget = NULL;

			// 个人理解：只在pivot之后处理相关行中会出现，如果新basic和 相关行 都有相同编号的Entry，就用一个指正指向，方便后续判断
            if ( _denseMap.exists( column ) )
                entryInTarget = _denseMap[column];

			//如果是处理新basic: 原basic行中的Entry都要 乘 scale变换系数（scale = ( -1.0 ) / 原basic的value值，只要乘完就是换算后的系数值）
			// 如果是处理Pivot后的相关行：scale是该行与新basic的系数值，假设相关行是新basic的3倍，则scale=3,而getValue()取得的是新basic与其他non-basic的倍数，
			// 两者相乘，得到的就是相关行相对于其他non-basic的倍数，（相当于就是将新basic带入到它此前曾出现过的其他等式中，计算系数）
            double newValue = current->getValue() * scale;

            // Statistics
            if ( numCalcs )
                ++( *numCalcs );

			// 个人理解：只有在pivot之后处理相关行会出现，即pivot之前，新basic还是non-basic在等式右边时，出现在了两个以上的等式右边，并且这两个等式有除了新baisc外其他相同的non-basic变量
            if ( entryInTarget )
            {
                if ( column != guaranteeIndex )
                    entryInTarget->setValue( entryInTarget->getValue() + newValue );	//不懂：为什么是+，不是直接替换
																					// 解答：这是处理在pivot之后处理相关行时，系数出现加减，原系数+-新系数
                else
                    entryInTarget->setValue( guaranteeValue );

                // The addition is another action
                if ( numCalcs )
                    ++( *numCalcs );

                if ( FloatUtils::isZero( entryInTarget->getValue() ) )
                {
                    _denseMap.erase( column );
                    eraseEntry( entryInTarget );
                }
            }
			// 个人理解：如果是处理新basic，一定会进入这个条件句子，生成新Entry
			// 如果是处理在pivot之后处理相关行，只有在本身没有时，才会引入新的，生成新Entry
            else
            {
                // This is a new entry for the target row
                if ( column != guaranteeIndex )
                    addEntry( target, column, newValue );
                else
                    addEntry( target, column, guaranteeValue );
            }
        }
    }

    void replaceNonBasicWithAnotherNonBasic(unsigned beReplaced, unsigned replace, double leakyValue){
        if ( !activeColumn( beReplaced ) )  // 将要被代换的，如果不存在，那就没有操作的必要，但将要去取代的，可以不存在
            return;

        _denseMap.clear();

        // 将 replace 出现过的行都记录下来，后面遍历beReplaced，取得出现过的行，如果两者有出现在同一行，就只是抹除beReplaced，
        // 然后在replace上做系数的加减，否则就要先抹除beReplaced,然后添加replace的系数
        Entry *replaceEntry = _columns[replace];
        while ( replaceEntry != NULL )
        {
            _denseMap[replaceEntry->getRow()] = replaceEntry;
            replaceEntry = replaceEntry->nextInColumn();
        }

        Entry *beReplacedEntry = _columns[beReplaced];
        Entry *current;

        unsigned row;
        while ( beReplacedEntry != NULL )
        {
            current = beReplacedEntry;
            beReplacedEntry = beReplacedEntry->nextInColumn();  // 记录列链表上下一个值，等待作为下一次while中判断条件

            row = current->getRow();
            Entry *entryInTarget = NULL;
            if ( _denseMap.exists( row ) )  // 如果要代替和被代替的变量出现在同一行上，改变系数、抹除在column链上的被代替的值
            {
                entryInTarget = _denseMap[row];
//                printf("~~~~~ The entryInTarget: row:%u, column: %u, value: %.10f\n", entryInTarget->getRow(), entryInTarget->getColumn(), entryInTarget->getValue());
                entryInTarget->setValue( entryInTarget->getValue() + leakyValue * current->getValue()); // 系数相加
//                printf("~~~~~ replaced and replace appear in the same row, set new coefficient value of %u: %.10f\n",entryInTarget->getRow(), entryInTarget->getValue());
            }
            else
            {
                addEntry(row, replace, leakyValue * current->getValue());  //添加系数，如果原来在x3f上的系数是1，这里要乘以0.5，表示换到x3b后是0.5，因为x3f=0.5x3b
            }
        }
        eraseColumn( beReplaced );  // 抹除所有beReplaced出现过的列
        _denseMap.clear();
    }

    void addColumnEraseSource( unsigned source, unsigned target )
    {
        if ( !activeColumn( source ) )
            return;

        _denseMap.clear();

        Entry *targetEntry = _columns[target];
        while ( targetEntry != NULL )
        {
            _denseMap[targetEntry->getRow()] = targetEntry;
            targetEntry = targetEntry->nextInColumn();
        }

        Entry *sourceEntry = _columns[source];
        Entry *current;

        unsigned row;
        while ( sourceEntry != NULL )
        {
            current = sourceEntry;
            sourceEntry = sourceEntry->nextInColumn();

            row = current->getRow();
            Entry *entryInTarget = NULL;
            if ( _denseMap.exists( row ) )
            {
                // Add from source column to target column
                entryInTarget = _denseMap[row];
                entryInTarget->setValue( entryInTarget->getValue() + current->getValue() );

                if ( FloatUtils::isZero( entryInTarget->getValue() ) )
                {
                    _denseMap.erase( row );
                    eraseEntry( entryInTarget );
                }
            }
            else
            {
                // There was no entry in the target column. "Steal" the entry.
                // No need to change row pointers, just the column pointers.

                if ( current->nextInColumn() )
                    current->nextInColumn()->setPrevInColumn( current->prevInColumn() );
                if ( current->prevInColumn() )
                    current->prevInColumn()->setNextInColumn( current->nextInColumn() );
                if ( _columns[source] == current )
                    _columns[source] = current->nextInColumn();

                // Adding is done at the head of the column, so shouldn't affect this loop
                current->setColumn( target );
                current->setNextInColumn( _columns[target] );
                current->setPrevInColumn( NULL );

                if ( _columns[target] )
                    _columns[target]->setPrevInColumn( current );

                _columns[target] = current;

                ++_columnSize[target];
                --_columnSize[source];
            }
        }

        eraseColumn( source );
        _denseMap.clear();
    }

    unsigned getRowSize( unsigned row ) const
    {
        return _rowSize.get( row );
    }

    unsigned getColumnSize( unsigned column ) const
    {
        return _columnSize.get( column );
    }

    const Entry *getRow( unsigned row ) const
    {
        return _rows[row];
    }

    const Entry *getColumn( unsigned column ) const
    {
        return _columns[column];
    }



    void printRow( unsigned row )
    {
        printf( "Printing row %u\n", row );
        printf( "\t%u = ", row );

        Entry *rowEntry = _rows[row];
        Entry *current;

        unsigned column;
        double weight;

        while ( rowEntry != NULL )
        {
            current = rowEntry;
            rowEntry = rowEntry->nextInRow();

            column = current->getColumn();
            weight = current->getValue();

            if ( !FloatUtils::isNegative( weight ) )
                printf( "+" );
            printf( "%lf * %u ", weight, column );
        }
        printf( "\n" );
    }

    void backupIntoMatrix( Tableau *other ) const
    {
        // Assume other has been initialized with correct sizes
        if ( other->_size != _size )
            throw Error( Error::COPY_INCOMPATIBLE_SPARSE_MATRICES );

        other->deleteAllEntries();
        other->_rowSize = _rowSize;
        other->_columnSize = _columnSize;

        for ( unsigned i = 0; i < _size; ++i )
        {
            Entry *entry = _rows[i];
            while ( entry )
            {
                Entry *newEntry = new Entry;

                newEntry->setRow( entry->getRow() );
                newEntry->setColumn( entry->getColumn() );
                newEntry->setValue( entry->getValue() );

                newEntry->setNextInRow( other->_rows[entry->getRow()] );
                newEntry->setNextInColumn( other->_columns[entry->getColumn()] );

                if ( other->_rows[entry->getRow()] )
                    other->_rows[entry->getRow()]->setPrevInRow( newEntry );
                if ( other->_columns[entry->getColumn()] )
                    other->_columns[entry->getColumn()]->setPrevInColumn( newEntry );

                other->_rows[entry->getRow()] = newEntry;
                other->_columns[entry->getColumn()] = newEntry;

                entry = entry->nextInRow();
            }
        }
    }

    void ensureNoZeros() const
    {
        for ( unsigned i = 0; i < _size; ++i )
            ensureNoZerosInRow( i );
    }

    void ensureNoZerosInRow( unsigned row ) const
    {
        if ( _rows[row] == NULL )
            return;

        Entry *entry = _rows[row];

        while ( entry != NULL )
        {
            if ( FloatUtils::isZero( entry->getValue() ) )
            {
                printf( "Error! Found a 0 in the matrix!\n" );
                exit( 1 );
            }

            entry = entry->nextInRow();
        }
    }

    unsigned countActiveColumns() const
    {
        unsigned result = 0;

        for ( unsigned i = 0; i < _columnSize.size(); ++i )
        {
            if ( _columnSize.get( i ) > 0 )
                result += 1;
        }

        return result;
    }



private:
    unsigned _size;
    Entry **_rows;	// 指针数组，初始化以后数组定长度，每个元素指向一个Entry指针，用于在AddEntry时记录双向链表的头，
					// 以便于构造双向链表，另一作用是在getRow方法中返回row双向链表的head,以便于遍历链表
    Entry **_columns;
    Vector<unsigned> _rowSize;
    Vector<unsigned> _columnSize;
    Map<unsigned, Entry *> _denseMap;

    void dumpDenseMap()
    {
        printf( "Dumping dense map (size = %u):\n", _denseMap.size() );
        Map<unsigned, Entry *>::iterator it = _denseMap.begin();
        while ( it != _denseMap.end() )
        {
            printf( "\t%u: %.5lf (address: 0x%p)\n", it->first, it->second->getValue(), it->second );
            ++it;
        }
    }
};

#endif // __Tableau_h__

//
// Local Variables:
// compile-command: "make -C . "
// tags-file-name: "./TAGS"
// c-basic-offset: 4
// End:
//
