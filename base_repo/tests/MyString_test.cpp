
/**************** Auto Generated File **********************/

#include "../includes/MyString.h"
#include <gtest/gtest.h>

const char kHelloString[] = "Hello, world!";

TEST(MyString, DefaultConstructor) {
  MyString s;
  EXPECT_STREQ(NULL, s.c_string());
  EXPECT_EQ(0u, s.Length());
}

TEST(MyString, Constructor_kHelloString) {
  MyString s(kHelloString);
  EXPECT_EQ(0, strcmp(s.c_string(), kHelloString));
  EXPECT_EQ(sizeof(kHelloString)/sizeof(kHelloString[0]) - 1,s.Length());
}

TEST(MyString, CopyConstructor) {
  MyString s1(kHelloString);
  MyString s2 = s1;
  EXPECT_EQ(0, strcmp(s2.c_string(), kHelloString));
}

TEST(MyString, Set) {
  MyString s;

  s.Set(kHelloString);
  EXPECT_EQ(0, strcmp(s.c_string(), kHelloString));

  s.Set(s.c_string());
  EXPECT_EQ(0, strcmp(s.c_string(), kHelloString));

  s.Set(NULL);
  EXPECT_STREQ(NULL, s.c_string());
}

