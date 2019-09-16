template<typename T> string tostr(T i)
{
	stringstream ss;
	ss << i;
	return ss.str();
}
