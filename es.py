class Solution(object):
    def merge(self, nums1, m, nums2, n):
        """
        :type nums1: List[int]
        :type m: int
        :type nums2: List[int]
        :type n: int
        :rtype: None Do not return anything, modify nums1 in-place instead.
        """
        i=0
        j=0
        l=[]
        while i<m or j<n:
            if i==m:
                l.append(nums2[j])
                j+=1
                continue
            if j==n:
                l.append(nums1[i])
                i+=1
                continue
            if (nums1[i]>=nums2[j]):
                l.append(nums2[j])
                j+=1
            else :
                l.append(nums1[i])
                i+=1
        return l
sol = Solution()
print(sol.merge([1,2,3,0,0,0],3,[2,5,6],3))