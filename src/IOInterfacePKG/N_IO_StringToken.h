

namespace Xyce {
namespace IO {

struct StringToken
{
  friend class Pack<StringToken>;

  StringToken(): lineNumber_(0), string_("")
  {}

  ~StringToken()
  {}

  size_t lineNumber_;
  std::string string_;
};

}
}