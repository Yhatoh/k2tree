#include<bits/stdc++.h>
using namespace std;

template <typename A, typename B>
string to_string(pair<A, B> p);
 
template <typename A, typename B, typename C>
string to_string(tuple<A, B, C> p);
 
template <typename A, typename B, typename C, typename D>
string to_string(tuple<A, B, C, D> p);

template<typename A, typename B, typename C, typename D, typename E>
string to_string(tuple<A, B, C, D, E> p);
 
inline string to_string(string s) {
  return '"' + s + '"';
}
inline string to_string(char s) {
  return '\'' + string(1, s) + '\'';
}
inline string to_string(const char* s) {
  return to_string((string) s);
}

inline string to_string(bool b) {
  return (b ? "true" : "false");
}
 
inline string to_string(vector<bool> v) {
  bool first = true;
  string res = "{";
  for (int i = 0; i < static_cast<int>(v.size()); i++) {
    if (!first) {
      res += ", ";
    }
    first = false;
    res += to_string(v[i]);
  }
  res += "}";
  return res;
}
 
template <size_t N>
string to_string(bitset<N> v) {
  string res = "";
  for (size_t i = 0; i < N; i++) {
    res += static_cast<char>('0' + v[i]);
  }
  return res;
}

template <typename A>
string to_string(queue<A> q){
  string res = "{";
  bool first = true;
  while(!q.empty()){
    if (!first)
      res += " <- ";
    res += to_string(q.front());
    first = false;
    q.pop();
  }
  res += "}";
  return res;
}

template <typename A>
string to_string(stack<A> q){
  string res = "[";
  stack<A> S;
  while(!q.empty()){
    S.push(q.top());
    q.pop();
  }

  bool first = true;
  while(!S.empty()){
    if (!first)
      res += " -> ";
    res += to_string(S.top());
    first = false;
    S.pop();
  }
  res += "}";
  return res;
}

template <class T, class Container,  class Compare>
string to_string(priority_queue<T, Container, Compare> q){
  string res = "{";
  bool first = true;
  while(!q.empty()){
    if (!first)
      res += " <- ";
    res += to_string(q.top());
    first = false;
    q.pop();
  }
  res += "}";
  return res;
}

template <typename A>
string to_string(A v) {
  bool first = true;
  string res = "{";
  for (const auto &x : v) {
    if (!first) {
      res += ", ";
    }
    first = false;
    res += to_string(x);
  }
  res += "}";
  return res;
}
 
template <typename A, typename B>
string to_string(pair<A, B> p) {
  return "(" + to_string(p.first) + ", " + to_string(p.second) + ")";
}
 
template <typename A, typename B, typename C>
string to_string(tuple<A, B, C> p) {
  return "(" + to_string(get<0>(p)) + ", " + to_string(get<1>(p)) + ", " + to_string(get<2>(p)) + ")";
}
 
template <typename A, typename B, typename C, typename D>
string to_string(tuple<A, B, C, D> p) {
  return "(" + to_string(get<0>(p)) + ", " + to_string(get<1>(p)) + ", " + to_string(get<2>(p)) + ", " + to_string(get<3>(p)) + ")";
}
 
template <typename A, typename B, typename C, typename D, typename E>
string to_string(tuple<A, B, C, D, E> p) {
  return "(" + to_string(get<0>(p)) + ", " + to_string(get<1>(p)) + ", " + to_string(get<2>(p)) + ", " + to_string(get<3>(p)) + ", " + to_string(get<4>(p)) + ")";
}

inline void debug_out() { cerr << endl; }
 
template <typename Head, typename... Tail>
void debug_out(Head H, Tail... T) {
  cerr << " " << to_string(H);
  debug_out(T...);
}
 
#define debug(...) cerr << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)
