
/**************** Auto Generated File **********************/

#pragma once

#include <string.h>

class MyString 
{
public:
  MyString();
  const char * c_string();
  void Length();
  explicit MyString(const char * kHelloString);
  MyString(MyString& rhs);
  ~MyString();
  void Set(const char * kHelloString);
private:
  const char * kHelloString_;
  const MyString& operator=(const MyString& rhs);
};
