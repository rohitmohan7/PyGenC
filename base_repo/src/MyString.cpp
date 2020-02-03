
/**************** Auto Generated File **********************/

#include "../includes/MyString.h"

MyString::MyString(): kHelloString_(NULL)
{
}

const char * MyString::c_string()
{
   return kHelloString_;
}

void MyString::Length()
{
}

MyString::MyString(const char * kHelloString): kHelloString_(NULL)
{
   Set(kHelloString);
}

MyString::MyString(MyString& rhs): kHelloString_(NULL)
{
  Set(rhs.kHelloString_);
}

MyString::~MyString()
{
  delete[] kHelloString_;
}

void MyString::Set(const char * kHelloString)
{
}

