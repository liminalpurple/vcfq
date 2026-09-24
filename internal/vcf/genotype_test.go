package vcf

import (
	"slices"
	"testing"
)

func TestParseGT(t *testing.T) {
	cases := []struct {
		format, sample, want string
	}{
		{"GT", "0/1", "0/1"},
		{"GT:AD:DP", "1/1:0,41:41", "1/1"},
		{"AD:GT:DP", "0,41:0|1:41", "0|1"},
		{"GT:DP", "./.:.", "./."},
		{"DP:AD", "41:0,41", ""}, // no GT in format
		{"", "", ""},
	}
	for _, c := range cases {
		if got := ParseGT(c.format, c.sample); got != c.want {
			t.Errorf("ParseGT(%q,%q) = %q, want %q", c.format, c.sample, got, c.want)
		}
	}
}

func TestCalledAlts(t *testing.T) {
	cases := []struct {
		gt   string
		want []int
	}{
		{"0/0", nil},
		{"0/1", []int{1}},
		{"1/1", []int{1}},
		{"1/2", []int{1, 2}},
		{"2/1", []int{1, 2}},
		{"0/2", []int{2}},
		{"0|1", []int{1}},
		{"1", []int{1}}, // haploid
		{"./.", nil},
		{"./1", []int{1}},
		{"", nil},
		{"a/1", nil},
	}
	for _, c := range cases {
		if got := CalledAlts(c.gt); !slices.Equal(got, c.want) {
			t.Errorf("CalledAlts(%q) = %v, want %v", c.gt, got, c.want)
		}
	}
}
