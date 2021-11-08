#include "stdafx.h"
#include "MASTER.h"
 #define ___3860
#include "GLOBAL.h"
#include "TASSERT.h"
#include "STRUTIL.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "STRLIST.h"
using tecplot::___1097; static ___372 StringListItemDestructor(void*      ___2098, ___90 ___494) { REQUIRE(VALID_REF(___2098)); REQUIRE(VALID_REF(*static_cast<char**>(___2098)) || *static_cast<char**>(___2098) == NULL); ___4278(___494); char** StringRef = static_cast<char**>(___2098); if (*StringRef != NULL) { ___1530(*StringRef, "string"); *StringRef = NULL; } ENSURE(*StringRef == NULL); return ___4226; } static ___372 StringListItemDuplicator(void*      ___3949, void*      ___3645, ___90 ___494) { REQUIRE(VALID_REF(___3949)); REQUIRE(VALID_REF(___3645)); REQUIRE(VALID_REF(*static_cast<char**>(___3645)) || *static_cast<char**>(___3645) == NULL); ___4278(___494); ___372 ___2040 = ___4226; char** TargetStringRef = static_cast<char**>(___3949); char** SourceStringRef = static_cast<char**>(___3645); if (*SourceStringRef != NULL) ___2040 = ((*TargetStringRef = ___1135(___1097(*SourceStringRef))) != NULL); else *TargetStringRef = NULL; ENSURE(VALID_REF(*TargetStringRef) || *TargetStringRef == NULL); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } ___372 ___3848(___3839 ___3819) { ___372 ___2067 = ArrayListIsValid(reinterpret_cast<___134>(___3819)); if (___2067) { ___2227 stringCount = ___101(reinterpret_cast<___134>(___3819));
 #if defined PERFORM_EXPENSIVE_STRLIST_TESTS
{ for (___2227 index = 0; index < stringCount; index++) { char* string = ___100(reinterpret_cast<___134>(___3819), index); if (string != NULL && !VALID_REF(string)) { ___2067 = ___1305; break; } } }
 #else
{ if (stringCount > 0) { char* string = ___100(reinterpret_cast<___134>(___3819), 0); if (string != NULL && !VALID_REF(string)) { ___2067 = ___1305; } } if (___2067 && stringCount > 1) { char* string = ___100(reinterpret_cast<___134>(___3819), stringCount - 1); if (string != NULL && !VALID_REF(string)) { ___2067 = ___1305; } } }
 #endif 
} ENSURE(VALID_BOOLEAN(___2067)); return ___2067; } void ___3824(___3839 ___3820) { REQUIRE(___3848(___3820)); ArrayListDeleteAllItems(reinterpret_cast<___134>(___3820), StringListItemDestructor, 0); ENSURE(___3848(___3820) && ___3826(___3820) == 0); } void ___3841(___3839 ___3820, ___2227     ___3854, ___2227     ___684) { REQUIRE(___3848(___3820)); REQUIRE(0 <= ___3854 && ___3854 <= ___3826(___3820) - 1); REQUIRE(1 <= ___684 && ___3854 + ___684 <= ___3826(___3820)); ArrayListDeleteItems(reinterpret_cast<___134>(___3820), ___3854, ___684, StringListItemDestructor, 0); ENSURE(___3848(___3820)); } void ___3840(___3839 ___3820, ___2227     ___3854) { REQUIRE(___3848(___3820)); REQUIRE(0 <= ___3854 && ___3854 <= ___3826(___3820) - 1); ArrayListDeleteItems(reinterpret_cast<___134>(___3820), ___3854, 1, StringListItemDestructor, 0); ENSURE(___3848(___3820)); } void ___3828(___3839* ___3820) { REQUIRE(VALID_REF(___3820)); REQUIRE(*___3820 == NULL || ___3848(*___3820)); if (*___3820 != NULL) ArrayListDealloc(reinterpret_cast<___134*>(___3820), StringListItemDestructor, 0); ENSURE(*___3820 == NULL); } ___3839 ___3821(void) { ___3839 ___3359 = reinterpret_cast<___3839>(ArrayListAlloc(0, ArrayListType_CharPtr, NULL, 0)); ENSURE(___3359 == NULL || ___3848(___3359)); return ___3359; } ___372 ___3823(___3839 ___3820, char const*   ___3813) { REQUIRE(___3848(___3820)); REQUIRE(___3813 == NULL || VALID_REF(___3813)); ___372 ___2040 = ___3843(___3820, ___3826(___3820), ___3813); ENSURE(___3848(___3820)); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } ___2227 ___3826(___3839 ___3820) { REQUIRE(___3848(___3820)); ___2227 ___3359 = ___101(reinterpret_cast<___134>(___3820)); ENSURE(___3359 >= 0); return ___3359; } char* ___3834(___3839 ___3820, ___2227     ___3854) { REQUIRE(___3848(___3820)); REQUIRE(0 <= ___3854 && ___3854 <= ___3826(___3820) - 1); char* ___3359; char const* StringRef = ___3835(___3820, ___3854); if (StringRef == NULL) ___3359 = NULL; else ___3359 = ___1135(___1097(StringRef)); ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; }
 #if !defined USE_MACROS_FOR_FUNCTIONS
char const* ___3836(___3839 ___3820, ___2227     ___3854) { REQUIRE(___3848(___3820)); REQUIRE(0 <= ___3854 && ___3854 <= ___3826(___3820) - 1); char const* ___3359 = ___3837(___3820, ___3854); ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; }
 #endif
___372 ___3843(___3839 ___3820, ___2227     ___3854, char const*   ___3813) { REQUIRE(___3848(___3820)); REQUIRE(___3854 >= 0); REQUIRE(___3813 == NULL || VALID_REF(___3813)); ___372       ___2040; ArrayListItem_u ItemCopy; if (___3813 != NULL) { ItemCopy.___474 = ___1135(___1097(___3813)); ___2040 = (ItemCopy.___474 != NULL); } else { ItemCopy.___474 = NULL; ___2040 = ___4226; } if (___2040) ___2040 = ArrayListSetItem(reinterpret_cast<___134>(___3820), ___3854, ItemCopy, StringListItemDestructor, 0); ENSURE(___3848(___3820)); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } ___372 ___3838(___3839 ___3820, ___2227     ___3854, char const*   ___3813) { REQUIRE(___3848(___3820)); REQUIRE(___3854 >= 0); REQUIRE(___3813 == NULL || VALID_REF(___3813)); ___372       ___2040; ArrayListItem_u ItemCopy; if (___3813 != NULL) { ItemCopy.___474 = ___1135(___1097(___3813)); ___2040 = (ItemCopy.___474 != NULL); } else { ItemCopy.___474 = NULL; ___2040 = ___4226; } if (___2040) ___2040 = ArrayListInsertItem(reinterpret_cast<___134>(___3820), ___3854, ItemCopy); ENSURE(___3848(___3820)); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } ___3839 ___3825(___3839 ___3820) { REQUIRE(___3848(___3820)); ___3839 ___3359 = reinterpret_cast<___3839>(ArrayListCopy(reinterpret_cast<___134>(___3820), StringListItemDuplicator, 0)); ENSURE(___3359 == NULL || (___3848(___3359) && ___3826(___3359) == ___3826(___3820))); return ___3359; } ___372 ___3822(___3839 ___3946, ___3839 ___3642) { REQUIRE(___3848(___3946)); REQUIRE(___3848(___3642)); ___3839 SourceCopy = ___3825(___3642); ___372 ___2040 = (SourceCopy != NULL); if (___2040) { ArrayListAppend(reinterpret_cast<___134>(___3946), reinterpret_cast<___134>(SourceCopy)); ArrayListDealloc(static_cast<___134*>(static_cast<void*>(&SourceCopy)), NULL, 0); } ENSURE(___3848(___3946)); ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } char* ___3847(___3839 ___3820) { REQUIRE(___3848(___3820)); size_t ___2224 = 0; ___2227 ___684 = ___3826(___3820); if (___684 >= 1) { ___2227 ___1926; for (___1926 = 0, ___2224 = strlen("\n") * (___684 - 1); ___1926 < ___684; ___1926++) { char* ___3813 = ___100(reinterpret_cast<___134>(___3820), ___1926); if (___3813 != NULL) ___2224 += strlen(___3813); } } char* ___3359 = ___23(___2224 + 1, char, "new line separated string"); if (___3359 != NULL) { ___2227 ___1926; for (___1926 = 0, strcpy(___3359, ""); ___1926 < ___684; ___1926++) { char* ___3813 = ___100(reinterpret_cast<___134>(___3820), ___1926); if (___1926 != 0) strcat(___3359, "\n"); if (___3813 != NULL) strcat(___3359, ___3813); } } ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; } ___3839 ___3831(char const* ___3813) { REQUIRE(VALID_REF(___3813)); ___3839 ___3359 = ___3821(); ___2227     ___3685; ___2227     EndIndex; for (___3685 = EndIndex = 0; ___3359 != NULL; EndIndex++) { if (___3813[EndIndex] == '\n' || ___3813[EndIndex] == '\0') { ___2227 ___2224 = EndIndex - ___3685; char*     SubString = ___23(___2224 + 1, char, "sub string"); if (SubString != NULL) { ___678(SubString, ___3813, ___3685, ___2224); ___3823(___3359, SubString); ___1530(SubString, "sub string"); if (___3813[EndIndex] != '\0') ___3685 = EndIndex + 1; else break; } else { ___3828(&___3359); ___3359 = NULL; break; } } } ENSURE(___3359 == NULL || ___3848(___3359)); return ___3359; } char** ___3845(___3839 ___3820) { REQUIRE(___3848(___3820)); char** ___3359 = static_cast<char**>(ArrayListToArray(reinterpret_cast<___134>(___3820), StringListItemDuplicator, 0)); ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; } ___3839 ___3829(char const** ___3816, ___2227    ___684) { REQUIRE((___684 == 0 && ___3816 == NULL) || (___684 >= 1 && VALID_REF(___3816))); ___3839 ___3359 = reinterpret_cast<___3839>(ArrayListFromArray(static_cast<void*>(___3816), ___684, ArrayListType_CharPtr, StringListItemDuplicator, 0)); ENSURE(___3359 == NULL || ___3848(___3359)); return ___3359; }
 #define ISJOINCHAR(c) ((c == ';') || (c == '+'))
static void SkipWhiteSpaceOrComma(char const** ___685) { REQUIRE(VALID_REF(___685) && VALID_REF(*___685)); while (___2082(**___685) || (**___685 == ',')) (*___685)++; } static ___372 GetNextSubString(char const** OriginalCPtr, char**       NextSubString) { REQUIRE(VALID_REF(OriginalCPtr) && (VALID_REF(*OriginalCPtr))); REQUIRE(VALID_REF(NextSubString)); ___372 ___2040 = ___4226; *NextSubString = NULL; char const* ___685 = *OriginalCPtr; SkipWhiteSpaceOrComma(&___685); char InsideDelimiter = '\0'; if (*___685 == '"'|| *___685 == '\'') { InsideDelimiter = *___685; ___685++; } char const* CStart = ___685; while (*___685 && ((InsideDelimiter && (*___685 != InsideDelimiter)) || (!InsideDelimiter && (*___685 != ',')       && !ISJOINCHAR(*___685)  && !___2082(*___685)))) { if (InsideDelimiter  && (*___685 == '\\')  && (*(___685 + 1) == InsideDelimiter)) ___685 += 2; else ___685++; } if (InsideDelimiter && (*___685 != InsideDelimiter)) ___2040 = ___1305; if (___2040 && CStart < ___685) { size_t StrLen = static_cast<size_t>(___685 - CStart); *NextSubString = ___23(StrLen + 1, char, "GetNextSubString: NextSubString"); if (*NextSubString) { char* NPtr = *NextSubString; while (CStart < ___685) { if ((*CStart == '\\') && (*(CStart + 1) == InsideDelimiter)) CStart++; *NPtr++ = *CStart++; } *NPtr = '\0'; } else ___2040 = ___1305; } if (___2040) { if (InsideDelimiter) ___685++; SkipWhiteSpaceOrComma(&___685); *OriginalCPtr = ___685; } ENSURE(VALID_BOOLEAN(___2040)); return ___2040; } ___3839 ___3830(char const* ___3813) { REQUIRE(VALID_REF(___3813)); SkipWhiteSpaceOrComma(&___3813); REQUIRE(!ISJOINCHAR(*___3813)); ___372 ___2040 = ___4226; ___3839 ___3359 = ___3821(); char const*   ___685   = ___3813; char* CurString = NULL; while (___2040 && *___685 != '\0') { char*     NextSubString = NULL; ___372 WantsToJoin   = ___1305; if (ISJOINCHAR(*___685)) { WantsToJoin = ___4226; ___685++; SkipWhiteSpaceOrComma(&___685); } ___2040 = GetNextSubString(&___685, &NextSubString); if (___2040) { if (WantsToJoin) ___3939(&CurString, '\n'); if (NextSubString != NULL && strlen(NextSubString) != 0) ___2040 = ___3941(&CurString, NextSubString, ___1305, ___1305); else if (CurString == NULL) CurString = ___1135(___1097("")); } if (NextSubString != NULL) ___1530(NextSubString, "StringListFromCompound: NextSubString"); if (___2040 && !ISJOINCHAR(*___685)) { ___3823(___3359, CurString); if (CurString != NULL) ___1530(CurString, "current string"); CurString = NULL; } } if (CurString != NULL) ___1530(CurString, "current string"); if (!___2040) ___3828(&___3359); ENSURE(___3359 == NULL || ___3848(___3359)); return ___3359; } char *___3846(___3839 ___3820, char          ___1818, char const*   ___475) { REQUIRE(___3848(___3820)); REQUIRE(___3826(___3820) >= 1); REQUIRE(ISJOINCHAR(___1818)); REQUIRE(VALID_REF(___475)); char* ___3359 = NULL; ___372 ___2040 = ___4226; ___2227 ___1926; ___2227 ___684; for (___1926 = 0, ___684 = ___3826(___3820), ___2040 = ___4226; ___1926 < ___684 && ___2040; ___1926++) { char* ___3813 = ___3834(___3820, ___1926); if (___3813 != NULL && strlen(___3813) != 0) { char*       CStart = NULL; char*       CEnd = NULL; char*       EscapedString = NULL; char const* EscChar = NULL; char*       StrChar = NULL; for (StrChar = ___3813; *StrChar != '\0'; StrChar++) { for (EscChar = ___475; *EscChar != '\0'; EscChar++) if (*StrChar == *EscChar) { ___2040 = ___3939(&EscapedString, '\\'); ___2040 = ___3939(&EscapedString, '\\'); break; } ___2040 = ___3939(&EscapedString, *StrChar); } CEnd = EscapedString; while (___2040 && CEnd && *CEnd != '\0') { int ___2223 = 0; CStart = CEnd; while (*CEnd != '\0' && *CEnd != '\n') { ___2223++; if (*CEnd =='"') ___2223++; CEnd++; } char* TString = ___23(___2223 + 4, char, "temp compound sub-string"); if (TString != NULL) { if (CStart == EscapedString) { if (___1926 != 0) ___2040 = ___3939(&___3359, ' '); } else { ___2040 = ___3939(&___3359, ___1818); } char* TStr = TString; *TStr++ ='"'; while (CStart && CStart != CEnd) { if (*CStart == '"') *TStr++ = '\\'; *TStr++ = *CStart++; } *TStr++ = '"'; *TStr = '\0'; ___3941(&___3359, TString, ___1305, ___1305); ___1530(TString, "___3846"); TString = NULL; if (*CEnd) CEnd++; } else { ___2040 = ___1305; } } if (EscapedString != NULL) ___1530(EscapedString, "escaped string"); } else { if (___1926 == 0) ___3941(&___3359, "\"\"", ___1305, ___1305); else ___3941(&___3359, " \"\"", ___1305, ___1305); } if (___3813 != NULL) ___1530(___3813, "string list ___2085"); } if (!___2040) { if (___3359 != NULL) { ___1530(___3359, "___3846"); ___3359 = NULL; } } ENSURE(___3359 == NULL || VALID_REF(___3359)); return ___3359; }