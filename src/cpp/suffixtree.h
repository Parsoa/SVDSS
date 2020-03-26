#ifndef _SUFFIX_TREE_HPP_INCLUDED_
#define _SUFFIX_TREE_HPP_INCLUDED_

#include <iostream>
#include <unordered_map>
#include <list>
#include <utility>
#include <memory>
#include <iterator>
#include <limits>
#include <vector>
#include <algorithm>

template <typename CharType = char, CharType end_token = '$'>
class SuffixTree {
    // Forward declarations of inner classes
    struct Node;

    typedef typename std::vector<CharType> string;
    typedef typename std::iterator_traits<typename string::iterator>::difference_type index_type;
    typedef CharType character;
    typedef std::tuple<Node*,index_type, index_type> ReferencePoint;


    // NESTED CLASSES DEFINITIONS

    // Node class:
    // Contains the suffix link to another Node
    // The Transitions "g(s,(k,i)) = s'" to children nodes
    // 
    // Note:
    // Transitions are stored in a hashtable indexed by the first substring
    // character. A given character can only have at most one Transition in a
    // node.
    
    // A Generalized Suffix Tree can contain more than one string at a time
    // Each string is labeled with an int. Thus each substring is related to
    // an appropriate reference string:
    // (ref string id, left ptr, right ptr)
    struct MappedSubstring {
        int ref_str;
        // A substring is a pair of integer (left ptr, right ptr)
        // To denote an empty substring, set right ptr < left ptr.
        index_type l;
        index_type r;
        MappedSubstring() : ref_str(0), l(0), r(0) {}
        MappedSubstring(int ref, index_type left, index_type right) :
          ref_str(ref),
          l(left),
          r(right)
          {}
        bool empty() const {
            return (this->l > this->r);
        }
    };

    struct Transition {
        MappedSubstring sub;
        Node *tgt;
        Transition() : sub(), tgt(nullptr) {}
        Transition(MappedSubstring s, Node *t) : sub(s), tgt(t) {}
        bool is_valid() {
            return (tgt != nullptr);
        }
    };

    struct Node {
        std::unordered_map<CharType, Transition> g;
        Node *suffix_link;
        virtual Transition find_alpha_transition(CharType alpha) {
            auto it = g.find(alpha);
            if (g.end() == it) {
                return Transition(MappedSubstring(0, 0, -1), nullptr);
            }
            return it->second;
        }
        
        Node() : suffix_link(nullptr) {}
        virtual ~Node() {}
        
        void dump_info() {
            for (auto t : g) {
                std::cout << "Transition for character: " << t.first << std::endl;
            }
        }
    };
    
    // Simple workaround for the specific sink node
    // This node must have a transition for every char of the input alphabet.
    // Instead of creating such transitions, we just make them up through
    // an override of `find_alpha_transition`
    struct SinkNode : public Node {
        virtual Transition find_alpha_transition(CharType alpha) override {
            return Transition(MappedSubstring(0, 0, 0), this->suffix_link);
        }
    };

    // Leaf nodes:
    // Leaves must contain an explicit reference to the suffix they represent
    // Some strings might have common suffixes, hence the map.
    // The suffix link **remains** UNIQUE nonetheless.
    struct Leaf : public Node {
        // TODO
    };

    // Base - A tree nested base class
    // This clase is here to hide implementation details
    // And to handle destruction properly.
    //
    // The processing (insertion, deletion of strings) is done by SuffixTree,
    // Base handles the cleanup.
    struct Base {
        SinkNode sink;
        Node root;
        Base() {
            root.suffix_link = &sink;
            sink.suffix_link = &root;
        }
        ~Base() {
            clean();
        }
        void clean() {
            std::list<Node*> del_list {&root};
            while (!del_list.empty()) {
                Node *current = del_list.front();
                del_list.pop_front();
                for (auto it : current->g) {
                    del_list.push_back(it.second.tgt);
                }
                if (&root != current) {
                    delete current;
                }
            }
        }
    };

    // "OUTER" CLASS MEMBERS

    Base tree;

    std::unordered_map<int, string> haystack;
    std::unordered_map<int, Node*> borderpath_map;
    int last_index;
    
    std::string to_string(string const & s, index_type b, index_type e) {
        std::string result;
        if (0 <= b && e < s.size()) {
            for (auto i = b; i <= e; ++i) {
                result.push_back(s[i]);
            }
        }
        return result;
    }
    
    std::string to_string(string const & s) {
        return to_string(s, 0, s.size()-1);
    }

    std::string to_string(MappedSubstring const & substr) {
        const auto& it = haystack.find(substr.ref_str);
        if (haystack.end() != it) {
            return to_string(it->second, substr.l, substr.r);
        }
    }

    // Given a Node n, a substring kp and a character t,
    // test_and_split must return if (n, kp) is the end point.
    // If not, and we are in an implicit state, a new state is created.
    bool test_and_split(Node *n, MappedSubstring kp, CharType t, const string& w, Node **r) {
        CharType tk = w[kp.l];
        index_type delta = kp.r - kp.l;
        if (0 <= delta) {
            Transition tk_trans = n->find_alpha_transition(tk);
            MappedSubstring kp_prime = tk_trans.sub;
            const auto& str_prime = haystack.find(kp_prime.ref_str);
            if (str_prime->second[kp_prime.l + delta + 1] == t) {
                *r = n;
                return true;
            } 
            *r = new Node();
            Transition new_t = tk_trans;
            new_t.sub.l += delta+1;
            (*r)->g.insert(std::pair<CharType, Transition>(
                str_prime->second[new_t.sub.l], new_t));
            tk_trans.sub.r = tk_trans.sub.l + delta;
            tk_trans.tgt = *r;
            n->g[tk] = tk_trans;
            return false;
                
        } else {
            // kp represents an empty substring
            Transition t_Transition = n->find_alpha_transition(t);
            *r = n;
            return (t_Transition.is_valid());
        }
    }

    // update performs the heart of an iteration:
    // It walks the border path from the active point to the end point
    // and adds the required Transitions brought by the insertion of
    // the string's i-th character.
    //
    // It returns the end point.
    ReferencePoint update(Node *n, MappedSubstring ki) {
        Node *oldr = &tree.root;
        Node *r = nullptr;
        bool is_endpoint = false;
        MappedSubstring ki1 = ki;
        const auto& ref_str_it = haystack.find(ki.ref_str);
        const string& w = ref_str_it->second;
        ReferencePoint sk(n, ki.ref_str, ki.l);
        ki1.r = ki.r-1;
        is_endpoint = test_and_split(n, ki1, w[ki.r], w, &r);
        while (!is_endpoint) {
            Leaf *r_prime = new Leaf();
            r->g.insert(std::make_pair(
              w[ki.r], Transition(MappedSubstring(
              ki.ref_str, ki.r, std::numeric_limits<index_type>::max()), r_prime)));
            if (&tree.root != oldr) {
                oldr->suffix_link = r;
            }
            oldr = r;
            sk = canonize(std::get<0>(sk)->suffix_link, ki1);
            ki1.l = ki.l = std::get<2>(sk);
            is_endpoint = test_and_split(std::get<0>(sk), ki1, w[ki.r], w, &r);
        }
        if (&tree.root != oldr) {
            oldr->suffix_link = std::get<0>(sk);
        }
        return sk;
    }

    // canonize - Get canonical pair
    // Given a Node and a substring,
    // This returns the canonical pair for this particular combination
    ReferencePoint canonize(Node *n, MappedSubstring kp) {
        if (kp.r < kp.l)
            return ReferencePoint(n, kp.ref_str, kp.l);
        const auto& kp_ref_str = haystack.find(kp.ref_str);
        index_type delta;
        Transition tk_trans = n->find_alpha_transition(kp_ref_str->second[kp.l]);
        while ((delta = tk_trans.sub.r - tk_trans.sub.l) <= kp.r - kp.l) {
            kp.l += 1 + delta;
            n = tk_trans.tgt;
            if (kp.l <= kp.r)
                tk_trans = n->find_alpha_transition(kp_ref_str->second[kp.l]);
        }
        return ReferencePoint(n, kp.ref_str, kp.l);
    }

    // get_starting_node - Find the starting node
    // @s[in]: The string to insert
    // @r[in/out]: The walk starting/ending point
    //
    // get_starting_node walks down the tree until s does not match anymore
    // character.
    // @r is updated through this process and contains the reference pair of the
    // diverging point between @s and the tree.
    // The result '(s,k)' of this function may then be used to resume the Ukkonen's
    // algorithm.
    index_type get_starting_node(const string& s, ReferencePoint *r) {
        auto k = std::get<2>(*r);
        auto s_len = s.size();
        bool s_runout = false;
        while (!s_runout) {
            Node *r_node = std::get<0>(*r);
            if (k >= s_len) {
                s_runout = true;
                break;
            }
            auto t = r_node->find_alpha_transition(s[k]);
            if (nullptr != t.tgt) {
                index_type i;
                const auto& ref_str = haystack.find(t.sub.ref_str);
                for (i=1; (i <= t.sub.r - t.sub.l); ++i) {
                    if (k+i >= s_len) {
                        s_runout = true;
                        break;
                    }
                    if (s[k+i] != ref_str->second[t.sub.l+i]) {
                        std::get<2>(*r) = k;
                        return k+i;
                    }
                }
                if (!s_runout) {
                   std::get<0>(*r) = t.tgt;
                   k += i;
                   std::get<2>(*r) = k;
                }
            } else {
                return k;
            }
        }
        std::get<2>(*r) = std::numeric_limits<index_type>::max();
        return std::numeric_limits<index_type>::max();
    }

    // deploy_suffixes - Deploy suffixes
    // @s[in]: The string to insert in the tree
    // @sindex[in]: The index id of @s
    //
    // deploy_suffixes performs the Ukkonen's algorithm to inser @s into the
    // tree.
    int deploy_suffixes(const string& s, int sindex) {
        ReferencePoint active_point(&tree.root, sindex, 0);
        auto i = get_starting_node(s, &active_point);
        if (std::numeric_limits<index_type>::max() == i) {
            return -1;
        }
        for (; i < s.size(); ++i) {
            MappedSubstring ki(sindex,std::get<2>(active_point), i);
            active_point = update(std::get<0>(active_point), ki);
            ki.l = std::get<2>(active_point);
            active_point = canonize(std::get<0>(active_point), ki);
        }
        return sindex;
    }

    void dump_node(Node *n, bool same_line, index_type padding, MappedSubstring orig) {
        index_type delta = 0;
        if (!same_line) {
            for (index_type i = 0; i < padding; ++i) {
                std::cout << "  ";
            }
        }
        if (!orig.empty()) {
            const auto& s = haystack.find(orig.ref_str);
            for (index_type i = orig.l; i <= orig.r && i < s->second.size(); ++i) {
                std::cout << s->second[i] << " ";
            }
            std::cout << "- ";
            delta = orig.r - orig.l + 2;
            if (orig.r == std::numeric_limits<index_type>::max()) {
                delta = s->second.size() - orig.l;
            }
        }
        same_line = true;
        for (auto t_it : n->g) {
            dump_node(t_it.second.tgt, same_line, padding + delta, t_it.second.sub);
            same_line = false;
        }
        if (same_line) {
            std::cout << "##" << std::endl;
        }
    }
    
    template <typename InputIterator>
    bool contain_end_token(InputIterator const & str_begin, InputIterator const & str_end) {
        return (std::find(str_begin, str_end, end_token) != str_end);
    }
    
    template <typename InputIterator, bool append_end_token = true>
    string make_string(InputIterator const & str_begin, InputIterator const & str_end) {
        if (contain_end_token(str_begin, str_end)) {
            throw std::invalid_argument("Input range contains the end token");
        }
        return make_string(str_begin, str_end, std::integral_constant<bool, append_end_token>());
    }
    
    template <typename InputIterator>
    string make_string(InputIterator const & str_begin, InputIterator const & str_end, std::true_type) {
        string s(str_begin, str_end);
        s.push_back(end_token);
        return s;
    }
    
    template <typename InputIterator>
    string make_string(InputIterator const & str_begin, InputIterator const & str_end, std::false_type) {
        index_type str_len = std::distance(str_begin, str_end);
        string s(str_begin, str_end);
        return s;
    }
public:
    SuffixTree() : last_index(0) {
    }
    
    template <typename InputIterator>
    int add_string(InputIterator const & str_begin, InputIterator const & str_end) {
        ++last_index;
        auto s = make_string(str_begin, str_end);
        haystack.emplace(last_index, std::move(s));
        const auto& s_from_map = haystack.find(last_index);
        if (0 > deploy_suffixes(s_from_map->second, last_index)) {
            haystack.erase(last_index--);
            return -1;
        }
        return last_index;
    }
    
    template <typename InputIterator>
    bool is_suffix(InputIterator const & str_begin, InputIterator const & str_end) {
        auto s = make_string(str_begin, str_end);
        ReferencePoint root_point(&tree.root, -1, 0);
        return (get_starting_node(s, &root_point) == std::numeric_limits<index_type>::max());
    }
    
    template <typename InputIterator>
    bool is_substring(InputIterator const & str_begin, InputIterator const & str_end) {
        auto s = make_string<false>(str_begin, str_end);
        ReferencePoint root_point(&tree.root, -1, 0);
        return (get_starting_node(s, &root_point) == std::numeric_limits<index_type>::max());
    }

    ~SuffixTree() {
    }

    void dump_tree() {
        dump_node(&tree.root, true, 0, MappedSubstring(0,0,-1));
    }
};

#endif // _SUFFIX_TREE_HPP_INCLUDED_
