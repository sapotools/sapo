/**
 * @file LinearSystemSet.h
 * Represent and manipulate a set of linear systems
 * It can be used to represent a symbolic union of polytopes
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEMSET_H_
#define LINEARSYSTEMSET_H_

#include "LinearSystem.h"

#define MINIMIZE_LS_SET_REPRESENTATION true


class LinearSystemSet {

private:
	using container = std::vector<LinearSystem*>;
	container set;	// set of linear systems

public:

	/**
	 * An iterator class for linear system set
	 */
	class iterator {
	
	private:
		container::iterator iit;

	public:
		using iterator_category = container::iterator::iterator_category;
		using difference_type = container::iterator::difference_type;
		using value_type = LinearSystem;
		using pointer = LinearSystem*;
		using reference = LinearSystem&;

		/**
		 *  Constructor
		 */
		inline iterator(const container::iterator& it): iit(it) {}

		/**
		 *  Copy constructor
		 */
		inline iterator(const iterator &orig): iit(orig.iit) {}

		/**
		 *  Swap constructor
		 */
		inline iterator(iterator&& orig) {
			std::swap(iit, orig.iit);
		}

		/**
		 *  Reference operator
		 *
		 * @return a reference to the pointed linear system.
		 */
		inline reference operator*() const
		{
			return *(*iit); 
		}

		/**
		 *  Pointer operator
		 *
		 * @return a pointer to the pointed linear system.
		 */
		inline pointer operator->() {
			return *iit;
		}

		/**
		 *  Prefix increment
		 *
		 * @return this iterator after incrementing it.
		 */
		inline iterator& operator++() {
			iit++; 
			return *this;
		}

		/**
		 *  Postfix increment
		 *
		 * @return this iterator before incrementing it.
		 */
		iterator operator++(int) {
			return iterator(iit++);
		}

		/**
		 *  Check whether two iterators point the same object.
		 * 
		 * @param[in] a is an iterator
		 * @param[in] b is an iterator
		 * @return true if and only if the two iterators point the same object.
		 */
		friend inline bool operator==(const iterator& a, const iterator& b) { 
			return a.iit == b.iit; 
		}

		/**
		 *  Check whether two iterators point the same object.
		 * 
		 * @param[in] a is an iterator
		 * @param[in] b is an iterator
		 * @return true if and only if the two iterators point different objects.
		 */
		friend inline bool operator!=(const iterator& a, const iterator& b) { 
			return a.iit != b.iit; 
		}
	};

	/**
	 * A constant iterator class for linear system set
	 */
	class const_iterator {
	private:
		container::const_iterator iit;

	public:
		using iterator_category = container::const_iterator::iterator_category;
		using difference_type = container::const_iterator::difference_type;
		using value_type = LinearSystem;
		using pointer = LinearSystem*;
		using reference = LinearSystem&;

		/**
		 *  Constructor
		 */
		inline const_iterator(const container::const_iterator& it): iit(it) {}

		/**
		 *  Copy constructor
		 */
		inline const_iterator(const const_iterator &orig): iit(orig.iit) {}

		/**
		 *  Swap constructor
		 */
		inline const_iterator(const_iterator&& orig) { 
			std::swap(iit, orig.iit);
		}

		/**
		 *  Reference operator
		 *
		 * @return a reference to the pointed linear system.
		 */
		inline reference operator*() const
		{
			return *(*iit); 
		}

		/**
		 *  Pointer operator
		 *
		 * @return a pointer to the pointed linear system.
		 */
		inline pointer operator->() {
			return *iit;
		}

		/**
		 *  Prefix increment
		 *
		 * @return this iterator after incrementing it.
		 */
		inline const_iterator& operator++() {
			iit++; 
			return *this;
		}

		/**
		 *  Postfix increment
		 *
		 * @return this iterator before incrementing it.
		 */
		const_iterator operator++(int) {
			return const_iterator(iit++);
		}

		/**
		 *  Check whether two iterators point the same object.
		 * 
		 * @param[in] a a constant iterator.
		 * @param[in] b a constant iterator.
		 * @return true if and only if the two iterators point the same object.
		 */
		friend inline bool operator==(const const_iterator& a, const const_iterator& b) { 
			return a.iit == b.iit; 
		}

		/**
		 *  Check whether two constant iterators point the same object.
		 * 
		 * @param[in] a is a constant iterator.
		 * @param[in] b is a constant iterator.
		 * @return true if and only if the two iterators point different objects.
		 */
		friend inline bool operator!=(const const_iterator& a, const const_iterator& b) { 
			return a.iit != b.iit; 
		} 
	};

	LinearSystemSet();
	LinearSystemSet(const LinearSystem& ls);
	LinearSystemSet(LinearSystem* ls);
	LinearSystemSet(const std::vector<LinearSystem*>& set);

	void add(LinearSystem *ls);
	void add(const LinearSystem& ls);
	void add(LinearSystem&& ls);

	LinearSystemSet& simplify();

	LinearSystemSet* get_a_finer_covering() const;

	// operations on set
	LinearSystemSet* getIntersectionWith(const LinearSystemSet *LSset) const;

	// in-place set operations
	LinearSystemSet& unionWith(LinearSystemSet *LSset);
	LinearSystemSet& boundedUnionWith(LinearSystemSet *LSset, const unsigned int bound);

	double boundingVol() const;

	/**
	 * Get the size of this set, i.e,
	 * the number of linear systems
	 *
	 * @returns number of linear systems in the set
	 */
	inline unsigned int size() const{ return this->set.size(); }

	/**
	 * Get the number of variables
	 * 
	 * @returns number of columns of linear systems in the set
	 */
	unsigned int dim() const;

	/**
	 * Get the set of linear systems
	 *
	 * @returns the current collection of linear systems
	 */
	inline const std::vector<LinearSystem*>& get_set() const { return set; }

	iterator begin() { iterator res(std::begin(set)); return res;}
	iterator end() { iterator res(std::end(set)); return res; }

	const_iterator cbegin() const { const_iterator res(std::begin(set)); return res;}
	const_iterator cend() const { const_iterator res(std::end(set)); return res; }

	bool isEmpty() const;

	/**
	 * Print the set of linear systems
	 */
	void print() const;

	/**
	 * Print the linear system set in Matlab format (for plotregion script)
	 * 
	 * @param[in] os is the output stream
	 * @param[in] color color of the polytope to plot
	 */
	void plotRegion(std::ostream& os=std::cout, const char color=' ') const;

	~LinearSystemSet();
};

std::ostream& operator<<(std::ostream& out, const LinearSystemSet& ls);

#endif /* LINEARSYSTEMSET_H_ */
