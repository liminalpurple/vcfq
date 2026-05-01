package vcf

import "testing"

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

func TestAltIndex(t *testing.T) {
	cases := []struct {
		gt   string
		want int
	}{
		{"0/0", 0},
		{"0/1", 1},
		{"1/1", 1},
		{"1/2", 2},
		{"2/3", 3},
		{"0|1", 1},
		{"./.", -1},
		{"", -1},
	}
	for _, c := range cases {
		if got := AltIndex(c.gt); got != c.want {
			t.Errorf("AltIndex(%q) = %d, want %d", c.gt, got, c.want)
		}
	}
}
